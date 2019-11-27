
#若无安装包，请从 <https://www.anaconda.com/download/#download>下载对应版本
bash .soft/Anaconda2-4.4.0-Linux-x86_64.sh
# 配置环境变量, 环境变量就是告诉系统起哪些目录寻找你输入的命令
# 远程登录，配置.bash_profile
nano ~/.bash_profile
export PATH=~/anaconda2/bin:${PATH}
source ~/.bash_profile
# 查看conda路径是否是自己安装的 
which conda 
# 给conda增加更多通道和国内镜像 （提高速度）
conda config --add channels conda-forge # Lowest priority
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels r # Optional
conda config --add channels defaults
conda config --add channels bioconda 
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
# Anocanda清华镜像
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ 
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ 
# # 清华通道, 最高优先级
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ 
conda config --set show_channel_urls yes
# # 显示已有的通道
conda config --get channels
# 安装需要的软件
conda install bwa samtools bedtools 
# 我们使用的conda是python3.6的运行环境，而RseQC只在python2.7环境中运行
# 所以需要新建一个py27环境，里面是python的2.7版本
# http://mp.weixin.qq.com/s/A4_j8ZbyprMr1TT_wgisQQ
conda create -n py27 python=2.7
source /anaconda3/bin/activate py27
conda install RseQC deeptools macs2
source /anaconda3/bin/deactivate py27
## 具体调用时，可以激活环境调用，也可以直接全路径调用
## 安装的环境在 ~/anaconda2/env/py27/下，可执行程序在~/anaconda2/env/py27/bin下
## 在 .bash_profile中新增下面这句话
export PATH=~/anaconda2/bin:${PATH}:~/anaconda2/env/py27/bin
# USCS tools安装
rsync -aP rsync://hgdownload.cse.ucsc.edu/genome/admin/exe/linux.x86_64/ ~/ucsc
## 在 .bash_profile中新增下面这句话
export PATH=~/anaconda2/bin:${PATH}:~/anaconda2/env/py27/bin:~/ucsc
# Picard安装
# http://broadinstitute.github.io/picard/
# 下载后解压，在解压目录下新建文件 picard.sh，内容如下
mem=$1
shift
java -Xmx${mem} -Djava.io.tmpdir=./ -jar $(dirname $(readlink -f $0))/picard.jar $@ VALIDATION_STRINGENCY=LENIENT
# 最后放到环境变量中
#
# Homer安装
# 从http://homer.ucsd.edu/homer/download.html下载
wget http://homer.ucsd.edu/homer/configureHomer.pl
# 安装homer2主程序
perl configureHomer.pl -install
# 安装小鼠注释文件 (比较大，需要较长时间)
perl configureHomer.pl -install mm10
# 基因组数据下载
wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/\
Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz -O mm10.fa.gz
gunzip mm10.fa.gz | sed 's/^/>chr/' | cut -f 1 -d ' ' >mm10.fa
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from mm10.chromInfo"  > mm10.chrom.size
# 构建索引，
bwa index -p bwa_mm10 mm10.fa
# 基因注释数据下载
wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz
gunzip -c Mus_musculus.GRCm38.90.gtf.gz | grep -v '^#' | sed 's/^/chr/' >mm10.gtf
# 获得bed12文件
# 这2个程序是UCSC下的软件，在刚才ucsc命令出现处已安装
gtfToGenePred -ignoreGroupsWithoutExons mm10.gtf mm10.gtf.50505050.pred
genePredToBed mm10.gtf.50505050.pred mm10.bed12
/bin/rm -f mm10.gtf.50505050.pred
# 获取转录起始位点
sed 's/"/\t/g' mm10.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {ensn=$10; symbol=$16; \
        if($7=="+") {start=$4-1; up=start; if(up<0) up=0; dw=start+1; \
        print $1,up, dw, ensn, symbol, $7;} else if($7=="-") {start=$5-1; up=start; \
        dw=start-1; if(dw<0) dw=0; print $1,dw,up,ensn,symbol,$7}}}' \
        | sort -k1,1 -k2,2n >mm10.tss.bed
# 基因bed文件
awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {split($9,a,";"); gene=a[1]; \
    split(gene, b, "\""); gene=b[2]; symbol=a[4]; split(symbol, c, "\""); \
    symbol=c[2]; print $1,$4-1,$5,gene,symbol,$7;}}' mm10.gtf \
    | sort -k1,1 -k2,2n >mm10.gene.bed
# 去掉长度小于500的区域，也可以用全部区域 
awk 'BEGIN{OFS=FS="\t" }{if($3-$2>500) print $0;}' mm10.gene.bed >mm10.gene_pc.bed
# 启动子位置
sed 's/"/\t/g' mm10.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {ensn=$10; symbol=$16; \
    if($7=="+") {start=$4-1; up=start-1000; if(up<0) up=0; dw=start+500; \
    print $1,up, dw, ensn, symbol, $7;} else \
    if($7=="-") {start=$5-1; up=start+1000; dw=start-500; \
    if(dw<0) dw=0; print $1,dw,up,ensn,symbol,$7}}}' | sort -k1,1 -k2,2n >mm10.promoter.bed
# 数据下载，完整的SRR列表在 SRR_list 文件中
# 一个数据有2个run，需要合并起来
# 如果只有一个run，也就是只有一个SRR号，则不需要合并
# 如果是双端测序，合并时需要考虑拆分后的 _1和_2
# --split-3 如果是双端测序则拆分
fastq-dump -v --split-3 --gzip SRR207071
fastq-dump -v --split-3 --gzip SRR207072
cat SRR207071.fastq.gz SRR207072.fastq.gz >MEF_CTCF.fq.gz
## 如果是双端，代码类似下方
# cat SRR207071_1.fastq.gz SRR207072_1.fastq.gz >MEF_CTCF_1.fq.gz
# cat SRR207071_2.fastq.gz SRR207072_2.fastq.gz >MEF_CTCF_2.fq.gz
# 删除下载的sra文件，节省空间
(cd ~/ncbi/public/sra/; /bin/rm -f SRR207071.sra SRR207072.sra)
# 指定使用4个线程
fastqc --extract -t 4 MEF_CTCF.fq.gz
# 去除接头和低质量碱基
# trimmomatic.sh SE -threads 10 -phred33 MEF_CTCF.fq.gz MEF_CTCF.trimmomatic.fq.gz \
  ILLUMINACLIP:MEF_CTCF.adaptor.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:17
# 新建一个文件夹，存储对应样品的比对结果
mkdir -p MEF_CTCF
# 开始比对，使用8个线程，-M是为了后续程序兼容，一般都使用这个参数
# bwa_mem10为索引文件，若不在当前目录需要指定全路径
bwa mem -M -t 8 ~/genome/bwa_mm10 MEF_CTCF.fq.gz | gzip >MEF_CTCF/MEF_CTCF.sam.gz
# 统计reads比对结果
# & 表示放入后台执行
bam_stat.py -i MEF_CTCF/MEF_CTCF.sam.gz >MEF_CTCF/MEF_CTCF.stat.xls &
# sam转bam
samtools view -F4 -q 1 -b MEF_CTCF/MEF_CTCF.sam.gz -o MEF_CTCF/MEF_CTCF.final.bam
# 比对reads按坐标排序，这是MarkDuplicates要求的输入文件格式
samtools sort -@ 8 -T /tmp/MEF_CTCF -o MEF_CTCF/MEF_CTCF.sortC.bam MEF_CTCF/MEF_CTCF.final.bam
# 去除PCR duplicates
picard.sh 15g MarkDuplicates INPUT=MEF_CTCF/MEF_CTCF.sortC.bam \
        OUTPUT=MEF_CTCF/MEF_CTCF.rmdup.bam METRICS_FILE=MEF_CTCF/MEF_CTCF.rmdup.log \
        CREATE_INDEX=true REMOVE_DUPLICATES=true
# 获取Duplication rate
tail -n 4 MEF_CTCF/MEF_CTCF.rmdup.log | cut -f 8 



# 基因组覆盖度计算
bedtools genomecov -ibam MEF_CTCF/MEF_CTCF.rmdup.bam -bga >MEF_CTCF/MEF_CTCF.genomecov.bdg
awk 'ARGIND==1{if($4>0) cover+=$3-$2;}ARGIND==2{genome+=$2;}END\
    {print "Genome regions have been sequenced compaared to full genome seq:", 
    cover/genome;}' MEF_CTCF/MEF_CTCF.genomecov.bdg ~/genome/mm10.chrom.sizes
awk 'ARGIND==1{if($4>0) {cover+$3-$2; depth=($3-$2)*$4;}}END\
    {print "Depth for covered regions:", depth/cover;}' \
    MEF_CTCF/MEF_CTCF.genomecov.bdg 
# reads分布计算
read_distribution.py -i MEF_CTCF/MEF_CTCF.final.bam -r ~/genome/mm10.gtf.bed12 \
    >MEF_CTCF/MEF_CTCF.read_distrib.xls
mkdir summary
plotFingerprint -b MEF_CTCF/MEF_CTCF.rmdup.bam --labels MEF_CTCF \
    --plotFile summary/MEF_CTCF.rmdup.fingerprint.pdf \
    --plotTitle "Fingerprints for All (removing PCR duplicates)" --skipZeros \
    --numberOfSamples 90000 --extendReads 200 --numberOfProcessors 2
plotFingerprint -b MEF_CTCF/MEF_CTCF.rmdup.bam MEF_H3K27ac/MEF_H3K27ac.rmdup.bam \
MEF_H3K4me1/MEF_H3K4me1.rmdup.bam MEF_H3K4me3/MEF_H3K4me3.rmdup.bam \
MEF_Input/MEF_Input.rmdup.bam MESC_CTCF/MESC_CTCF.rmdup.bam \
MESC_H3K27ac/MESC_H3K27ac.rmdup.bam MESC_H3K4me1/MESC_H3K4me1.rmdup.bam \
MESC_H3K4me3/MESC_H3K4me3.rmdup.bam MESC_Input/MESC_Input.rmdup.bam --labels \
MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac \
MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotFile summary/All.rmdup.fingerprint.pdf \
--plotTitle "Fingerprints for All (removing PCR duplicates)" --skipZeros \
    --numberOfSamples 90000 --extendReads 200 --numberOfProcessors 10
bamCoverage -b MEF_CTCF/MEF_CTCF.rmdup.bam -o MEF_CTCF/MEF_CTCF.rmdup.coverage.bw \
    --binSize 50 --blackListFileName ~/genome/mm10.blacklist.bed \
    --ignoreForNormalization chrX chrM --normalizeUsing RPGC --smoothLength 100 \
    --numberOfProcessors 20 --effectiveGenomeSize 1870000000 \
    --extendReads 200 --centerReads
# IP-input
bamCompare -b1 MEF_CTCF/MEF_CTCF.rmdup.bam -b2 MEF_Input/MEF_Input.rmdup.bam \
    --outFileFormat bedgraph -o MEF_CTCF/MEF_CTCF.rmdup.substract.bdg \
    --blackListFileName ~/genome/mm10.blacklist.bed --operation subtract \
    --extendReads 200 --ignoreForNormalization chrX chrM --normalizeUsing RPKM \
    --centerReads --numberOfProcessors 20 
awk 'BEGIN{OFS=FS="\t"}{if($4<0) $4=0; print $0}' MEF_CTCF/MEF_CTCF.rmdup.substract.bdg \
    >MEF_CTCF/MEF_CTCF.rmdup.substract.bdg2
bedGraphToBigWig MEF_CTCF/MEF_CTCF.rmdup.substract.bdg2 ~/genome/mm10.chrom.sizes \
    MEF_CTCF/MEF_CTCF.rmdup.substract.bw
/bin/rm -f MEF_CTCF/MEF_CTCF.rmdup.substract.bdg*
touch All.gene.profile
mkdir -p summary
computeMatrix scale-regions --regionsFileName ~/genome/mm10.gene_pc.bed \
    --scoreFileName MEF_CTCF/MEF_CTCF.rmdup.substract.bw \
MEF_H3K27ac/MEF_H3K27ac.rmdup.substract.bw MEF_H3K4me1/MEF_H3K4me1.rmdup.substract.bw \
MEF_H3K4me3/MEF_H3K4me3.rmdup.substract.bw MESC_CTCF/MESC_CTCF.rmdup.substract.bw \
MESC_H3K27ac/MESC_H3K27ac.rmdup.substract.bw \
MESC_H3K4me1/MESC_H3K4me1.rmdup.substract.bw MESC_H3K4me3/MESC_H3K4me3.rmdup.substract.bw \
    --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
    --skipZeros -o summary/All.gene.profile.rmdup.substract.gz \
    --blackListFileName ~/genome/mm10.blacklist.bed --binSize 50 --missingDataAsZero \
    --numberOfProcessors 20
plotProfile --matrixFile summary/All.gene.profile.rmdup.substract.gz \
    --outFileName summary/All.gene.profile.rmdup.substract.meta.pdf --perGroup \
    --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac \
MESC_H3K4me1 MESC_H3K4me3 \
    --plotTitle "Meta-gene profile for all samples (control and PCR duplicates removed)" \
    --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotProfile --matrixFile summary/All.gene.profile.rmdup.substract.gz \
    --outFileName summary/All.gene.profile.rmdup.substract.meta.kmeans.pdf \
    --perGroup --kmeans 2 \
    --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac \
    MESC_H3K4me1 MESC_H3K4me3 \
    --plotTitle "Meta-gene profile for all samples (control and PCR duplicates removed)" \
    --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotHeatmap --matrixFile summary/All.gene.profile.rmdup.substract.gz \
    --outFileName summary/All.gene.profile.rmdup.substract.heatmap.kmeans.pdf \
    --kmeans 4 --colorMap Blues --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 \
    --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac \
    MESC_H3K4me1 MESC_H3K4me3 \
    --plotTitle "All.gene.profile (control and PCR duplicates removed)" \
    --heatmapWidth 6 --heatmapHeight 20
computeMatrix scale-regions --regionsFileName ~/genome/mm10.gene_pc.bed --scoreFileName MEF_CTCF/MEF_CTCF.rmdup.coverage.bw MEF_H3K27ac/MEF_H3K27ac.rmdup.coverage.bw MEF_H3K4me1/MEF_H3K4me1.rmdup.coverage.bw MEF_H3K4me3/MEF_H3K4me3.rmdup.coverage.bw MEF_Input/MEF_Input.rmdup.coverage.bw MESC_CTCF/MESC_CTCF.rmdup.coverage.bw MESC_H3K27ac/MESC_H3K27ac.rmdup.coverage.bw MESC_H3K4me1/MESC_H3K4me1.rmdup.coverage.bw MESC_H3K4me3/MESC_H3K4me3.rmdup.coverage.bw MESC_Input/MESC_Input.rmdup.coverage.bw --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o summary/All.gene.profile.rmdup.coverage.gz --blackListFileName ~/genome/mm10.blacklist.bed --binSize 50 --missingDataAsZero --numberOfProcessors 20
plotProfile --matrixFile summary/All.gene.profile.rmdup.coverage.gz --outFileName summary/All.gene.profile.rmdup.coverage.meta.pdf --perGroup --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "Meta-gene profile for all samples (PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotProfile --matrixFile summary/All.gene.profile.rmdup.coverage.gz --outFileName summary/All.gene.profile.rmdup.coverage.meta.kmeans.pdf --perGroup --kmeans 2 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "Meta-gene profile for all samples (PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotHeatmap --matrixFile summary/All.gene.profile.rmdup.coverage.gz --outFileName summary/All.gene.profile.rmdup.coverage.heatmap.kmeans.pdf --kmeans 4 --colorMap Blues --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "All.gene.profile (PCR duplicates removed)" --heatmapWidth 6 --heatmapHeight 32
# TSS profile
mkdir -p summary
computeMatrix reference-point --regionsFileName ~/genome/mm10.tss.bed --scoreFileName MEF_CTCF/MEF_CTCF.rmdup.substract.bw MEF_H3K27ac/MEF_H3K27ac.rmdup.substract.bw MEF_H3K4me1/MEF_H3K4me1.rmdup.substract.bw MEF_H3K4me3/MEF_H3K4me3.rmdup.substract.bw MESC_CTCF/MESC_CTCF.rmdup.substract.bw MESC_H3K27ac/MESC_H3K27ac.rmdup.substract.bw MESC_H3K4me1/MESC_H3K4me1.rmdup.substract.bw MESC_H3K4me3/MESC_H3K4me3.rmdup.substract.bw --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -o summary/All.tss.profile.rmdup.substract.gz --blackListFileName ~/genome/mm10.blacklist.bed --binSize 50 --missingDataAsZero --numberOfProcessors 20 
plotProfile --matrixFile summary/All.tss.profile.rmdup.substract.gz --outFileName summary/All.tss.profile.rmdup.substract.meta.pdf --perGroup --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 --plotTitle "Meta-gene profile for all samples (control and PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotProfile --matrixFile summary/All.tss.profile.rmdup.substract.gz --outFileName summary/All.tss.profile.rmdup.substract.meta.kmeans.pdf --perGroup --kmeans 2 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 --plotTitle "Meta-gene profile for all samples (control and PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotHeatmap --matrixFile summary/All.tss.profile.rmdup.substract.gz --outFileName summary/All.tss.profile.rmdup.substract.heatmap.kmeans.pdf --kmeans 4 --colorMap Blues --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 --plotTitle "TSS profile (control and PCR duplicates removed)" --heatmapWidth 6 --heatmapHeight 32
computeMatrix reference-point --regionsFileName ~/genome/mm10.tss.bed --scoreFileName MEF_CTCF/MEF_CTCF.rmdup.coverage.bw MEF_H3K27ac/MEF_H3K27ac.rmdup.coverage.bw MEF_H3K4me1/MEF_H3K4me1.rmdup.coverage.bw MEF_H3K4me3/MEF_H3K4me3.rmdup.coverage.bw MEF_Input/MEF_Input.rmdup.coverage.bw MESC_CTCF/MESC_CTCF.rmdup.coverage.bw MESC_H3K27ac/MESC_H3K27ac.rmdup.coverage.bw MESC_H3K4me1/MESC_H3K4me1.rmdup.coverage.bw MESC_H3K4me3/MESC_H3K4me3.rmdup.coverage.bw MESC_Input/MESC_Input.rmdup.coverage.bw --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -o summary/All.tss.profile.rmdup.coverage.gz --blackListFileName ~/genome/mm10.blacklist.bed --binSize 50 --missingDataAsZero --numberOfProcessors 20 
plotProfile --matrixFile summary/All.tss.profile.rmdup.coverage.gz --outFileName summary/All.tss.profile.rmdup.coverage.meta.pdf --perGroup --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "Meta-gene profile for all samples (PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotProfile --matrixFile summary/All.tss.profile.rmdup.coverage.gz --outFileName summary/All.tss.profile.rmdup.coverage.meta.kmeans.pdf --perGroup --kmeans 2 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "Meta-gene profile for all samples (PCR duplicates removed)" --yAxisLabel "Ratio" --plotWidth 15 --plotHeight 9
plotHeatmap --matrixFile summary/All.tss.profile.rmdup.coverage.gz --outFileName summary/All.tss.profile.rmdup.coverage.heatmap.kmeans.pdf --kmeans 4 --colorMap Blues --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --samplesLabel MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MEF_Input MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 MESC_Input --plotTitle "TSS profile (PCR duplicates removed)" --heatmapWidth 6 --heatmapHeight 32
multiBigwigSummary bins --bwfiles MEF_CTCF/MEF_CTCF.rmdup.coverage.bw MEF_H3K27ac/MEF_H3K27ac.rmdup.coverage.bw MEF_H3K4me1/MEF_H3K4me1.rmdup.coverage.bw MEF_H3K4me3/MEF_H3K4me3.rmdup.coverage.bw MESC_CTCF/MESC_CTCF.rmdup.coverage.bw MESC_H3K27ac/MESC_H3K27ac.rmdup.coverage.bw MESC_H3K4me1/MESC_H3K4me1.rmdup.coverage.bw MESC_H3K4me3/MESC_H3K4me3.rmdup.coverage.bw --labels MEF_CTCF MEF_H3K27ac MEF_H3K4me1 MEF_H3K4me3 MESC_CTCF MESC_H3K27ac MESC_H3K4me1 MESC_H3K4me3 --outFileName summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.npz --outRawCounts summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.tab --binSize 1000 --blackListFileName ~/genome/mm10.blacklist.bed --numberOfProcessors 20
plotCorrelation --corData summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.npz --plotFile summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.correlation.pdf --corMethod spearman --whatToPlot heatmap --skipZeros --plotTitle "Sample correlation (PCR duplicates removed)" --removeOutliers --plotNumbers --outFileCorMatrix summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.correlation.xls
sed "1 s/'//g" summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.tab | awk 'BEGIN{OFS=FS="\t"}{$1=$1"_"$2; print $0;}' | cut -f 2,3 --complement >summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.xls
s-plot pca -f summary/YSX_train.rmdup.coverage.multiBigwigSummary_bin.xls -T 50000 -L TRUE -N 3 -t "PCA after PCR duplicates and removed"
## sharp peak
macs2 callpeak -t MEF_CTCF/MEF_CTCF.final.bam -g 1870000000 --outdir MEF_CTCF -n MEF_CTCF.peak_call -B --SPMR -q 0.01 --call-summits --fix-bimodal -f AUTO --seed 11521 -c MEF_Input/MEF_Input.final.bam 
cut -f 1-3 MEF_CTCF/MEF_CTCF.peak_call_peaks.narrowPeak >MEF_CTCF/MEF_CTCF.peak.bed
## Broad peak
macs2 callpeak -t MESC_H3K4me1/MESC_H3K4me1.final.bam -g 1870000000 --outdir MESC_H3K4me1 -n MESC_H3K4me1.peak_call -B --SPMR --fix-bimodal --extsize 500 --broad --broad-cutoff 0.1 -f AUTO --seed 11521 -c MESC_Input/MESC_Input.final.bam 
cut -f 1-3 MESC_H3K4me1/MESC_H3K4me1.peak_call_peaks.narrowPeak >MESC_H3K4me1/MESC_H3K4me1.peak.bed
intersectBed -a MEF_CTCF/MEF_CTCF.peak_call_summits.bed -b ~/genome/mm10.promoter.bed -wa -wb >MEF_CTCF/MEF_CTCF.peak_gene_promoter.xls
cat MEF_CTCF/MEF_CTCF.peak_gene_promoter.xls | cut -f 9 | cut -d '.' -f 1 | sed 's/$/\tMEF_CTCF/'>MEF_CTCF/MEF_CTCF.peak_gene_promoter.list
bedtools closest -D a -fu -t all -a MEF_CTCF/MEF_CTCF.peak_call_summits.bed -b ~/genome/mm10.gene.bed | sed '1 i\chr\tstart_0_based\tend\tpeak_summit\tpeak_score\tchr\tstart_0_based\tend\tENSM\tSymbol\tStrand\tDistance_negvalue_upstream' >MEF_CTCF/MEF_CTCF.peak_gene_closest.xls
tail -n +2 MEF_CTCF/MEF_CTCF.peak_gene_closest.xls | cut -f 9 | cut -d '.' -f 1 | sed 's/$/\tMEF_CTCF/'>MEF_CTCF/MEF_CTCF.peak_gene_closest.list
# 在给定区域鉴定motif
# findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]
# -p: 多线程
# -size: 200 (中心200) , -100,50 上游100下游50，given给定区域
# -bg: homer会自动选取，也可以自己指定
findMotifsGenome.pl MEF_CTCF/MEF_CTCF.peak_call_summits.bed mm10 MEF_CTCF_Motif -p 12 -size 50 
bedtools getfasta -fi ~/genome/mm10.fa -bed MESC_CTCF/MESC_CTCF.peak.bed >MESC_CTCF/MESC_CTCF.peak.fa
bedtools shuffle -i MESC_CTCF/MESC_CTCF.peak.bed -g ~/genome/mm10.chrom.sizes -excl MESC_CTCF/MESC_CTCF.peak.bed >MESC_CTCF/MESC_CTCF.shuffle.bed
findMotifs.pl MESC_H3K4me1/MESC_H3K4me1.peak.fa fasta OutputDir -p 8 -fastaBg background.fa
# 在给定区域搜寻motif
findMotifsGenome.pl MESC_H3K4me1/MESC_H3K4me1.peak.bed mm10 -nomotif -find /data/chip/soft/data/knownTFs/all/all.motifs
#The result is output_dir/*.find.tmp
transferFindMotifsOutputToBed.py -i input/*.find.tmp -b input.bed -t DNA -H 0 >input.knownTfs.bed
# 在给定序列搜寻motif
bedtools getfasta -fi ~/genome/mm10.fa -bed MESC_H3K4me1/MESC_H3K4me1.peak.bed >MESC_H3K4me1/MESC_H3K4me1.peak.fa
findMotifs.pl MESC_H3K4me1/MESC_H3K4me1.peak.fa fasta OutputDir -p 8 -find /data/chip/soft/data/knownTFs/all/all.motifs >H3K4me1.knownTfs
transferFindMotifsOutputToBed.py -i H3K4me1.knownTfs -b input.bed -t DNA >H3K4me1.knownTfs.bed



