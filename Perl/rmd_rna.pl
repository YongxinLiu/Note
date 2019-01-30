#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

sub usage {
    die(
        qq!
Usage:    rmd_rna.pl -o output_dir -d design.txt -c group_compare.txt -v group_venn.txt
Function: write RNA-Seq report in Rbookdown
Command:  -i inpute file name (Must)
          -o output file directory (Must), defaul current directory
          -d design.txt
          -c group_compare.txt
          -v group_venn.txt
          -V report.txt
          -r report.txt
          -h header line number, default 0
          -p pvalue
          -q qvalue
          -S output SNP result, default FALSE
          -C output cross species result, default FLASE
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.3
Update:   2017/6/2
Notes:    1.0 output html report; 1.1 add GO and KEEG of clusterProfiler; 1.2 Add SNP and cross species; 1.3 debug GO/KEEG
\n!
    )
}

my %opts;
getopts( 'i:o:d:e:l:c:v:h:n:r:s:p:q:S:C:V:', \%opts );
$opts{h}=0 unless defined($opts{h});
$opts{o}="./" unless defined($opts{o});
$opts{d}="doc/design.txt" unless defined($opts{d});
$opts{c}="lsls" unless defined($opts{c});
$opts{v}="doc/group_venn.txt" unless defined($opts{v});
$opts{r}="doc/summary.txt" unless defined($opts{r});
$opts{e}="TRUE" unless defined($opts{e});
$opts{S}="FALSE" unless defined($opts{S});
$opts{C}="FALSE" unless defined($opts{C});
$opts{p}=0.05 unless defined($opts{p});
$opts{q}=0.2 unless defined($opts{q});

&usage unless (exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

# Prepare relative files
print "Clean report enviroment and prepare relative files\n";
`rm -f *.Rmd`;
`rm -fr _bookdown_files`;
`rm -fr $opts{V}`;
`cp -f -r /mnt/bai/yongxin/ref/RNA/rmd/* $opts{o}`;
`sed -i 's/html/$opts{V}/g' _bookdown.yml`;

my %list;
open LIST,"<$opts{r}";
while (<LIST>) {
	chomp;
	my @tmp=split/\t/;
	$list{$tmp[0]}=$tmp[1];
}
close LIST;

# Set project default parameter
$list{title}='转录组测序数据分析(RNA-Seq)' unless defined($list{title});
$list{partner}=$list{client} unless defined($list{partner});
$list{analyst}="Dr. Yong-Xin Liu, Bailab, IGDB, CAS" unless defined($list{analyst});
$list{contact}="yxliu\@genetics.ac.cn 010-64808722" unless defined($list{contact});
$list{project}="170901-1" unless defined($list{project});
$list{period}="2017-09-01~ 2017-12-31" unless defined($list{period});
$list{website}="http://bailab.genetics.ac.cn/" unless defined($list{website});
$list{wechat}="[植物微生物组](http://mp.weixin.qq.com/s/QkgNlzK_rpauKzSs2fUd7A)" unless defined($list{wechat});
$list{logo}="Bailab, SKLPG/CEPAMS, IGDB, CAS {-}" unless defined($list{logo});
$list{logo_file}="bailab" unless defined($list{logo_file});

open OUTPUT,">$opts{o}index.Rmd";
print OUTPUT qq!--- 
title: "$list{title}"
author:
- Partner实验数据：$list{partner}
- Analyst生信分析：$list{analyst}
- Contact联系方式：$list{contact}
- Project项目编号：$list{project}
- Period项目周期：$list{period}
- Website官方网站：$list{website}
- Wechat公众号：$list{wechat}
date: '`r Sys.Date()`'
documentclass: article
bibliography: [rna.bib]
link-citations: yes
biblio-style: apalike
---

```{r setup, include=FALSE}
library(knitr)
output <- opts_knit\$get("rmarkdown.pandoc.to")
html = FALSE
latex = FALSE
opts_chunk\$set(echo = FALSE, out.width="100%", fig.align="center", fig.show="hold", warning=FALSE, message=FALSE)
if (output=="html") {
	html = TRUE
}
if (output=="latex") {
	opts_chunk\$set(out.width="95%", out.height='0.7\\textheight', out.extra='keepaspectratio', fig.pos='H')
	latex = TRUE
}
knitr::opts_chunk\$set(cache=TRUE, autodep=TRUE)
mtime <- function(files){
  lapply(Sys.glob(files), function(x) file.info(x)\$mtime)
}
set.seed(718)
```

```{asis, echo=html}
# Bailab, SKLPG/CEPAMS, IGDB, CAS {-}
```

```{r cover, eval=html, out.width="99%"}
figs_1 = paste0("figure/banner", ".png")
knitr::include_graphics(figs_1)
```
!;
close OUTPUT;



open OUTPUT,">$opts{o}01-aim.Rmd";
print OUTPUT qq!
# 课题目的 {#project_aim}

$list{aim}

!;
close OUTPUT;



open OUTPUT,">$opts{o}02-design.Rmd";
print OUTPUT qq!
# 课题设计 {#project_design}

$list{design}，样品准备如 Table \\\@ref(tab:design) 所示。[design.txt]($opts{d})

```{r design}
table_design <- read.table("$opts{d}", sep="\\t", header=T)
knitr::kable(table_design, caption="样品详细及分组信息总结。", booktabs=TRUE)
```

!;
close OUTPUT;



open OUTPUT,">$opts{o}03-scheme.Rmd";
print OUTPUT qq!
# 课题方案 {#project_scheme}

**1. 测序reads数量和质量评估；Quality control of sequencing reads**

Table: (\\#tab:seq-quality-explanatioan-ch) 测序质量评估结果解读方法

-----------------------------------------------------------------------------------
评估内容                   结果解释 (图例中会标记对应评估内容为PASS、WARN和FAIL, 具体处理方式详见下面中英文解释)
-------------------------  --------------------------------------------------------------------------------
Per base quality           测序reads从5'到3'的碱基的质量值 (Q)。该值越大越代表对应碱基测序准确度越高。假设p为一个碱基测序错误的概率，则Q=-10 * log10(p). 质量值为10时，对应碱基出错的概率为10%；质量值为20时，对应碱基出错的概率为1%。通常来讲，3'端的碱基质量会低于5'端；另外5'端最初几个碱基也会出现较大的质量值波动。我们在后期处理时，会去除低质量的碱基以保证分析结果的准确性。

Adaptor content            判断测序reads中是否残留接头序列。存在接头序列和不存在接头序列都是合理的，取决于测序数据下机后是否进行了接头去除和去除的是否完整。若在分析时检测到接头序列存在，我们会首先去除接头，然后进行后续分析，以保证分析结果的准确性。

Per sequence GC content    测序reads的GC含量。正常的测序reads的GC含量符合正态分布模式 (形如图中蓝色的倒钟形线)。若在平滑的曲线上存在一个尖峰表示测序样品存在特定的序列污染如混入了引物二聚体。若GC含量分布曲线比较平坦则代表可能存在不同物种的序列污染。当这一指标异常时，可能导致后期的序列比对或拼接存在问题，需要引起注意。

Per base sequence content  测序reads的碱基偏好性。正常的测序结果中一个序列不同的碱基没有偏好性，图中的线应平行。Bisulfite测序中存在甲基化的C到T的转变，会导致这一评估结果异常。我们人工核验无误后，可以忽略软件对这一检测结果的评价。
-----------------------------------------------------------------------------------


Table: (\\#tab:seq-quality-explanatioan-en) Explanation of quality control by fastqc.

-----------------------------------------------------------------------------------
Analysis                   Explanation
-------------------------  --------------------------------------------------------------------------------
Per base quality           The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

Per sequence GC content    Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example),  which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.

Adaptor content            Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This doesn't indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.

Per base sequence content  In a random library you would expect that there would be little to no difference between the different bases of a sequence run,  so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome,  but in any case they should not be hugely imbalanced from each other.
-----------------------------------------------------------------------------------


(ref:scheme-read-fastqc) 测序Reads质量评估。HiSeqX产出Clean reads左端(A)和右端(B)各150 bp数据质量评估，选取测序reads碱基质量分布判断建库或测序有无异常; Reads左端(C)和右端(D)各碱基A/G/C/T含量分布，用于分析各样品测序物种的GC偏好性，有无接头及引物污染，接头去除干净与否以及有效数据比例评估；双端序列重复情况(E/F)，常用于评估物种的基因组倍性情况，重复序列，PCR扩增重复度等情况。  Quality control of raw reads. (A/B) Per base sequence quality; (C/D) Per base sequence content; (E/F) Sequence duplication levels. [\@andrews2010fastqc \@bolger2014trimmomatic]

```{r scheme-read-fastqc, fig.cap="(ref:scheme-read-fastqc)", out.width="49%"}
figs_1 = paste0("figure/", c("1.1quality", "1.2quality", "1.3content", "1.4content", "1.5duplication", "1.6duplication"),".png")
knitr::include_graphics(figs_1)
```

**2. 样品提取及过滤各步骤统计；Statistics of reds filter processes**

(ref:scheme-read-summary) 统计样品处理过程及样品可用数据量，如raw read, clean reads, mapped, gene, intergenic, intron, exon, UTR。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B). 柱状图展示各样品的数据量分布，最小值也大于2万，大部分在12万左右，完全符合实验设计要求；(C) 可用数据的长度分布，可以看到本实验扩增子长度范围集中在360-390 bp，主峰位于370-380 bp间。RNA-Seq常用PE150质量非常高，在有参考基因组的分析基本可以跳过质控步骤；无参考基因组的转录组则必须进行极严格的质量控制。  Statistics of reads filter processes in libraries and data size of samples. (A) Bar plot showing reads count of each library in read filter process; (B) Bar plot showing reads counts of each sample; (C) Length distribution of amplicons [\@bolger2014trimmomatic].

```{r scheme-read-summary, fig.cap="(ref:scheme-read-summary)"}
knitr::include_graphics("figure/fig2.summary.png")
```

**3. 样品基于参考基因组可定位比例分布。Statistics of mapping processes on reference genome**

(ref:scheme-read-map) 统计各样品与参考基因组比对结果中各类分布比例。PEunqiue, PEmultiple代表双端序列在基因组上可定位于参考基因组且位置符合预期(concordantly)，并且分别在基因组上有唯一位置(unique)和多个位置(multiple)；PEdiscordantly代表双端序列可定位于参考基因组，但位置不符合预期，如双端方向与预期相反、位置相距太远等，可能由于染色体变异或融合基因转录而成；SEunique和SEmultiple表示双端序列只有一端可比对到参考基因组上，分别为位置唯一或多个位置；Unmapped代表没有比对到参考基因组上，原因可能为参考基因组不完整或样品中微生物污染或其它污染导致。 [\@pertea2016transcript \@li2009sequence \@kim2015hisat \@dobin2013star].

```{r scheme-read-map, fig.cap="(ref:scheme-read-map)"}
knitr::include_graphics("figure/3mapping_stat.png")
```

**4. 样品和组相关分析 **

(ref:scheme-read-cor) 计算各样品和组间的Pearson相关系数，显示组间相似或差异关系，样品间重复性好坏的评估[\@anders2014htseq \@mckenna2010genome ]。

```{r scheme-read-cor, fig.cap="(ref:scheme-read-cor)", out.width="99%"}
figs_1 = paste0("figure/", c("4.1cor_groups", "4.2cor_samples"),".png")
knitr::include_graphics(figs_1)
```

**5. 差异表达基因分析；Differentially expressed genes**

(ref:scheme-sample-gene) KO1基因型存在一些丰度显著上调或下调的OTU (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示KO与WT相比OTU的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调，图中数字代表显著差异OTU数量，形状代表OTU的门水平物种注释；（B）热图展示KO与WT显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低；可以看到我们找到的差异OTU在每组样品中重复非常好，同时也发现了在beta diversity分析中发现的KO1中存在的两个异常样品应该为KO1.7, KO1.8, 需检查实验材料准备了取材步骤有无问题？或补弃样品重复（C）曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。  KO1 are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs of KO1 compared WT; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs. These results show Actinobacterial has more enriched OTUs [\@bai2015functional, \@edwards2015structure, \@zgadzaj2016root].

```{r scheme-sample-gene, fig.cap="(ref:scheme-sample-gene)"}
knitr::include_graphics("figure/fig8.otu.png")
```

**6. 差异基因的GO和KEGG富含分析 **

(ref:scheme-read-go) 计算各组差异基因的富含通路，包括基因本体论(Gene Ontology, GO)的生物学过程(Biological Process, BP)，分子功能(Molecular Function, MF)，细胞组分(Celluar Conponent, CC)，以及基因和基因组京都百科全书(Kyoto Encyclopedia of Genes and Genomes, KEGG)注释的基因通路。[\@yu2012clusterprofiler]

```{r scheme-read-go, fig.cap="(ref:scheme-read-go)", out.width="99%"}
figs_1 = paste0("figure/5.", c("1BP", "2MF", "3CC", "4KK"),".png")
knitr::include_graphics(figs_1)
```

(ref:scheme-read-gonet) 富含GO term详细分析。(A)GO条目(terms)包括相关基因图(centplot)，展示核心GO条目具体包含那些功能基因；(B)富含GO条目相互关系网络图(enrichMap)，展示GO条目间的相互作用关系；(C) GO条目层级关系图(GOgraph)，展示富集的GO条目在注释层级中的位置，其中颜色代表其pvalue值，越深代表越显著富集。[\@yu2012clusterprofiler \@wickham2016ggplot2]

```{r scheme-read-gonet, fig.cap="(ref:scheme-read-gonet)", out.width="99%"}
figs_1 = paste0("figure/5.", c("5cnetplot", "6enrichMap", "7GOgraph"),".png")
knitr::include_graphics(figs_1)
```

**7. 组间差异表达基因比较；Compare differentially expressed genes among groups**

(ref:scheme-sample-overlap) 比较组间差异基因。(A) 饼形图展示各种差异基因分类比例。中间数字为所有显著差异基因的数目。可以看到KO1与KO2样式相似，OE1与OE2样式相似。且上调OTU较多为Actinobacteria，而下调OTU绝大多数为Proteobacteria。(B) 维恩图展示各基因型差异OTUs间的共有和特有数量。图中所显各基因型组间重复间大部分OTUs共有；而且还发现KO和OE还存在一些相似变化样式的OTUs。  Taxonomy, common and unique OTUs in each group. (A) Pie charts show phyla of bacterial OTUs identified as either enriched or depleted in each genotype compared with WT. The number of OTUs in each category is noted inside each donut. (B) Venn diagrams show common and unique OTUs in each group [\@robinson2010edger \@love2014moderated \@chen2011venndiagram \@lebeis2015salicylic].

```{r scheme-sample-overlap, fig.cap="(ref:scheme-sample-overlap)"}
knitr::include_graphics("figure/fig9.overlap.png")
```

**9. 目标基因表达结果展示；Showing one gene expression**

(ref:scheme-gene-exp) 展示目标单个基因的表达情况。(A) 箱线图展示单基因在各组中的表达水平，点代表单个样品；(B) 基因结构热图展示基因在组内的平均表达水平；(C) 基因结构热图展示基因在每个样品中的表达水平。(A) Boxplot showing expression of single gene in different groups, point represents each sample; (B) Heatmap of gene structure showing the average expression value in group; (C) Heatmap of gene structure showing the gene expression value in each sample [\@frazee2015ballgown \@pertea2015stringtie].

```{r scheme-gene-exp, fig.cap="(ref:scheme-gene-exp)"}
figs_1 = paste0("figure/bg_", c("boxplot", "structure_groups", "structure_samples"),".png")
knitr::include_graphics(figs_1)
```

**9. 其它数据分析过程中发现的有意思的点，商讨后，有意义的深入分析；Other points and ideas for further discussion and analysis **

如亲本间差异表达基因分析[\@nodine2012maternal]；

!;
close OUTPUT;


## 读取样品列表文件
open DATABASE,"<$opts{d}";
<DATABASE>;
	while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @sample,$tmp[0];
}
close DATABASE;

open OUTPUT,">$opts{o}04-a-sequenceQuality.Rmd";
print OUTPUT qq!
# 测序质量总结 {#sequencing_quality_summary}

## 测序质量评估 {#sub-sequence-qc}

说明：RNA-Seq通常采用HiSeqX PE150测序，兼顾读长和价格，性价比超值。

!;

foreach $sample (@sample) {
#`unzip -o seq/${sample}.clean.R1_fastqc.zip -d seq/`;
#`unzip -o seq/${sample}.clean.R2_fastqc.zip -d seq/`;
print OUTPUT qq!
### 样品${sample}质量评估
(ref:quality-fastqc-${sample}) 测序Reads质量评估：样品${sample}。Clean reads左端(A)和右端(B)数据质量评估；clean reads左端(C)和右端(D)序列重复情况分布。Quality control of clean reads [HTML_1](seq/${sample}_1_fastqc.html)  [HTML_2](seq/${sample}_2_fastqc.html)

```{r quality-fastqc-${sample}, fig.cap="(ref:quality-fastqc-${sample})", out.width="49%"}
figs_1 = paste0("seq/${sample}_", c("1_fastqc/Images/per_base_quality", "2_fastqc/Images/per_base_quality", "1_fastqc/Images/duplication_levels", "2_fastqc/Images/duplication_levels"),".png")
knitr::include_graphics(figs_1)
```
!;
}

print OUTPUT qq!
## 样品比对参考基因组统计


```{r table-map}
table_fpkm <- read.table("result/statmap.txt", sep="\t", header=T,row.names=1)
table_fpkm = table_fpkm/rowSums(table_fpkm)*100
knitr::kable(table_fpkm, caption="样品比对参考基因组统计表[TXT](result/statmap.txt)", booktabs=TRUE)
```

(ref:quality-map) 测序数据与参考基因组比对统计图。Barplot show mapping percentage. [PDF](result/statmap.txt.fillBars.pdf)

```{r quality-map, fig.cap="(ref:quality-map)", out.width="99%"}
figs_1 = paste0("result/statmap.txt.fillBars", ".png")
knitr::include_graphics(figs_1)
```
!;
close OUTPUT;


open OUTPUT,">$opts{o}04-d-diversity.Rmd";
### 样品与组间表达值 {#result-exp}
#
#表 \\\@ref(tab:table-fpkm)展示了所有样品Top30基因表达FPKM值
#
#```{r table-fpkm}
#table_fpkm <- read.table("result/gene_fpkm_all.txt", sep="\t", header=T)
#table_fpkm = head(table_fpkm,n=30)
#knitr::kable(table_fpkm, caption="样品Top30基因表达FPKM值", booktabs=TRUE)
#```
#
#[所有样品基因FPKM值](result/gene_fpkm.txt)  [基因>5FPKM的样品及组值](result/gene_fpkm_all.anno)  [所有样品转录本FPKM值](result/transcript_fpkm.txt)  
#[所有样品基因测序reads统计](result/htseq_anno.txt)  
#
#
### 基因和转录本表达分布箱线图 {#result-exp-box}
#
#(ref:gene-box) 各样品基因表达箱线图，基因表达FPKM值采用以2为底的对数进行标准化。Boxplot showing distribution of gene expression log2(FPKM) . [PDF](result/gene_fpkm_boxplot.pdf)
#
#```{r gene-box, fig.cap="(ref:gene-box)", out.width="99%"}
#knitr::include_graphics("result/gene_fpkm_boxplot.png")
#```
#
#(ref:transcript-box) 各样品转录本表达箱线图，转录本表达FPKM值采用以2为底的对数进行标准化。Boxplot showing distribution of transcript expression log2(FPKM) . [PDF](result/transcript_fpkm_boxplot.pdf)
#
#```{r transcript-box, fig.cap="(ref:transcript-box)", out.width="99%"}
#knitr::include_graphics("result/transcript_fpkm_boxplot.png")
#```

print OUTPUT qq!
# 样品与组表达值及相关分析 {#result-diversity}


## 样品或组间Pearson相关分析 {#result-diversity-cor}

(ref:group-cor) 基于各样品组基因平均RKM值计算皮尔森相关系数。Pearson correlation of all groups. [PDF](result/heat_cor_groups.pdf)

```{r group-cor, fig.cap="(ref:group-cor)", out.width="99%"}
knitr::include_graphics("result/heat_cor_groups.png")
```

(ref:sample-cor) 基于各样品基因RPM值计算皮尔森相关系数。Pearson correlation of all samples. [PDF](result/heat_cor_samples.pdf)

```{r sample-cor, fig.cap="(ref:sample-cor)", out.width="99%"}
knitr::include_graphics("result/heat_cor_samples.png")
```
!;
close OUTPUT;



open OUTPUT,">$opts{o}04-f-gene.Rmd";
print OUTPUT qq!
#	差异基因分析 {#result-gene}

## 各组差异表达基因数量概述 {#result-gene-sum}

样品组间显著差异基因数量(P < $opts{p}, FDR < $opts{q})如下表所示。[数据Summary](result/gene_sum.txt) [详细列表Detail](result/gene.txt)

```{r gene-sum}
table_gene <- read.table("result/gene_sum.txt", sep="\t", header=T)
knitr::kable(table_gene, caption="各组差异表达基因数量汇总", booktabs=TRUE)
```


## 差异基因分析 {#result-gene-de}

!;

open DATABASE,"<$opts{c}";
while (<DATABASE>) {
	chomp;
	push @group,$_;
}
$i=0;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	$i++;
	`rm -f result/gene_$tmp[0]vs$tmp[1]_all.xls`;
	`cat result/gene_$tmp[0]vs$tmp[1]_*.xls > result/gene_$tmp[0]vs$tmp[1]_all.xls`;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

(ref:gene-$i) $tmp[0]vs$tmp[1]组间存在一些表达显著上调或下调的基因 (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示$tmp[0]与$tmp[1]相比基因表达的差异分布，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调；(B) 热图展示$tmp[0]与$tmp[1]显著差异基因在每个样品中相对丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示基因延染色体的差异分布情况，对于观察某些成簇表达变化的基因，x轴为基因在染色体上排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，大小为相对丰度，颜色和形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of number of d`ferentially expressed genes, fold change and expression level; (B) Heatmap showing significantlly  differentially expressed genes; (C) Manhattan plot showing differentially expressed genes along the chromosome.
[Volcano plot PDF](result/vol_gene_$tmp[0]vs$tmp[1].pdf)  [Heatmap PDF](result/heat_gene_$tmp[0]vs$tmp[1]_sig.pdf)  [Manhattan plot PDF](result/man_gene_$tmp[0]vs$tmp[1].pdf)

```{r gene-$i, fig.cap="(ref:gene-$i)", out.width="99%"}
figs_2 = paste0("result/", c("vol_gene_$tmp[0]vs$tmp[1]", "heat_gene_$tmp[0]vs$tmp[1]_sig", "man_gene_$tmp[0]vs$tmp[1]"),".png")
knitr::include_graphics(figs_2)
```
$tmp[0]与$tmp[1]相比显著差异基因信息。[Enriched](result/gene_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted](result/gene_$tmp[0]vs$tmp[1]_depleted.xls)  [No significant](result/gene_$tmp[0]vs$tmp[1]_nosig.xls)  [All](result/gene_$tmp[0]vs$tmp[1]_all.xls)

```{r tab-gene-$tmp[0]vs$tmp[1]}
e = read.table("result/gene_$tmp[0]vs$tmp[1]_enriched.txt", sep="\t", header=T)
e = head(e[order(-e\$A_mean),],n=10)
d = read.table("result/gene_$tmp[0]vs$tmp[1]_depleted.txt", sep="\t", header=T)
d = head(d[order(-d\$B_mean),],n=10)
m = rbind(e,d)
m=m[,c("ID","A_mean","B_mean","logFC","PValue","level")]
knitr::kable(m, row.names = F, caption="样品$tmp[0]与$tmp[1]相比显著上调或下调表达Top10基因(按表达丰度降序排列)；Significantlly different genes.", booktabs=TRUE)
```

$tmp[0]vs$tmp[1]组间存在一些表达显著上调或下调的基因功能富含分析(超几何分布检验Pvalue < $opts{p}, qvalue < $opts{q})。主要的分类有GO生物学过程、分子功能分、细胞组分和KEGG。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (Pvalue < $opts{p}, qvalue < $opts{q}, hypergeometric distribution test). Mainly included Enriched GO terms involved in biological process, molecular function, celluar component and KEEG pathway.

!;

#(ref:GO-$i) $tmp[0]vs$tmp[1]组间存在一些表达显著上调或下调的基因功能富含分析。 (P & FDR < 0.01, 超几何分布)。(A) GO生物学过程分类富集结果；(B) GO分子功能分类富集结果；(C) GO细胞组分分类富集结果。
#Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.01, hypergeometric distribution test). (A) Enriched GO terms involved in biological process; (B) Enriched GO terms involved in molecular function; (C) Enriched GO terms involved in celluar component; (D) Enriched KEEG pathway.
#[Enriched GO BP](result/gene_$tmp[0]vs$tmp[1]_enriched_BP.pdf)  [Depleted GO BP](result/gene_$tmp[0]vs$tmp[1]_depleted_BP.pdf)
#[Enriched GO MF](result/gene_$tmp[0]vs$tmp[1]_enriched_MF.pdf)  [Depleted GO MF](result/gene_$tmp[0]vs$tmp[1]_depleted_MF.pdf)
#[Enriched GO CC](result/gene_$tmp[0]vs$tmp[1]_enriched_CC.pdf)  [Depleted GO CC](result/gene_$tmp[0]vs$tmp[1]_depleted_CC.pdf)
#[Enriched KEEG](result/gene_$tmp[0]vs$tmp[1]_enriched_KK.pdf)  [Depleted KEEG](result/gene_$tmp[0]vs$tmp[1]_depleted_KK.pdf)
#[Enriched GO BP_enrichMap](result/gene_$tmp[0]vs$tmp[1]_enriched_BP_enrichMap.pdf)  [Depleted GO BP_enrichMap](result/gene_$tmp[0]vs$tmp[1]_depleted_BP_enrichMap.pdf)
#[Enriched GO BP_cnetplot](result/gene_$tmp[0]vs$tmp[1]_enriched_BP_cnetplot.pdf)  [Depleted GO BP_cnetplot](result/gene_$tmp[0]vs$tmp[1]_depleted_BP_cnetplot.pdf)
#[Enriched GO BP_plotGOgraph](result/gene_$tmp[0]vs$tmp[1]_enriched_BP_plotGOgraph.pdf)  [Depleted GO BP_plotGOgraph](result/gene_$tmp[0]vs$tmp[1]_depleted_BP_plotGOgraph.pdf)
#
#```{r GO-$i, fig.cap="(ref:GO-$i)", out.width="49%"}
#figs_2 = paste0("result/gene_$tmp[0]vs$tmp[1]_", c("enriched_BP", "depleted_BP","enriched_MF", "depleted_MF","enriched_CC", "depleted_CC","enriched_KK", "depleted_KK"),".png")
#knitr::include_graphics(figs_2)
#```

# 下滑线也会引起表头名字对出错
@tax_en=qw#enriched_BP enriched_MF enriched_CC enriched_KK depleted_BP depleted_MF depleted_CC depleted_KK#;
foreach $j (@tax_en) {
if (-e "result/gene_$tmp[0]vs$tmp[1]_$j.pdf") {

print OUTPUT qq!


(ref:GO-${i}-$j)[$j](result/gene_$tmp[0]vs$tmp[1]_$j.pdf)

```{r GO-${i}-$j, fig.cap="(ref:GO-${i}-$j)", out.width="99%"}
figs_2 = paste0("result/gene_$tmp[0]vs$tmp[1]_", c("$j"),".png")
knitr::include_graphics(figs_2)
```
!;
}
}

}

## 读group venn文件
open DATABASE,"<$opts{v}";
while (<DATABASE>) {
	chomp;
	push @venn,$_;
}
close DATABASE;

if (@venn>=1) {

print OUTPUT qq!
## 比较共有差异基因 {#result-gene-venn}

!;

$i=0;
foreach (@venn) {
	chomp;
	$i++;
	my @tmp=split/\t/; # sampleA and sampleB
	my $venn_list2;
	foreach $tmp (@tmp) {
		$venn_list2.=$tmp;
	}
	$tmp[2]="C" unless defined($tmp[2]);
	$tmp[3]="D" unless defined($tmp[3]);
	$tmp[4]="E" unless defined($tmp[4]);
	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";
	$id=$venn_list2;

print OUTPUT qq!

### $venn_list

(ref:gene-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result/gene.txt.venn$venn_list.pdf)  [List XLS](result/gene.txt.venn$venn_list2.xls)  [Detail XLSX](result/gene.txt.venn$venn_list2.xls.xls)  


```{r gene-venn-$i, fig.cap="(ref:gene-venn-$i)", out.width="99%"}
figs_2 = paste0("result/gene.txt.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```


!;


#print "result/gene.txt.venn${id}.xls.list\n";
$file = "result/gene.txt.venn${id}.xls.list";
if (-e $file) {
open DATABASE,"<result/gene.txt.venn${id}.xls.list";
while (<DATABASE>) {
	chomp;
	$id=$_;

	print $id,"\n";
print OUTPUT qq!

### ${id}

(ref:GO-${id}) ${id}组间存在一些表达显著上调或下调的基因功能富含分析。 (P & FDR < 0.01, 超几何分布)。(A) GO生物学过程分类富集结果；(B) GO分子功能分类富集结果；(C) GO细胞组分分类富集结果。
Compared with $tmp[1], $tmp[0] are and depleted for certain genes. (P & FDR < 0.01, hypergeometric distribution test). (A) GO terms involved in biological process; (B) GO terms involved in molecular function; (C) GO terms involved in celluar component; (D) KEEG pathway.
[GO BP](result/gene_${id}_BP.pdf)
[GO MF](result/gene_${id}_MF.pdf)
[GO CC](result/gene_${id}_CC.pdf)
[KEEG](result/gene_${id}_KK.pdf)
[GO BP_enrichMap](result/gene_${id}_BP_enrichMap.pdf)
[GO BP_cnetplot](result/gene_${id}_BP_cnetplot.pdf)
[GO BP_plotGOgraph](result/gene_${id}_BP_plotGOgraph.pdf)

```{r GO-${id}, fig.cap="(ref:GO-${id})", out.width="49%"}
figs_2 = paste0("result/gene_${id}_", c("BP", "MF", "CC", "KK"),".png")
knitr::include_graphics(figs_2)
```

!;
}
}
}
}












if ($opts{S} eq "TRUE") {
open OUTPUT,">$opts{o}04-g-snp.Rmd";
print OUTPUT qq!
# SNP分型确定等位基因差异表达 {#result-snp}

实验采用拟南芥Col和Cvi品种杂交，花粉来自Cvi，Cvi的单体型图谱来自拟南芥1000 genome；  
方法主要采用STAR2 mapping RNA序列于参考基因组上，再用samtools统计每个位点上的各种碱基数据；
基于Cvi的SNP位点信息，筛选每个样品中这些位置的reads count，进行统计可区分位点的整体和单基因表达情况。


## 样品中来自亲本双方的SNP位点统计 {#result-snp-sum}

```{r table-snp-count}
snp_sum_count <- read.table("result/snp_sum_count.txt", sep="\t", header=T, row.names = 1)
per=(snp_sum_count)/rowSums(snp_sum_count)*100
all=cbind(snp_sum_count,per)
knitr::kable(all, caption="SNP位点来自不同生态型中高质量reads数量及比例汇总", booktabs=TRUE)
```

SNP位点来自不同生态型中高质量reads数量[TXT](result/snp_sum_count.txt)


(ref:table-snp-bar) 各样品中基因分型比较柱状图 [PDF](result/snp_sum_count.txt.fillBars.pdf)

```{r table-snp-bar, fig.cap="(ref:table-snp-bar)", out.width="99%"}
figs_2 = paste0("result/", c("snp_sum_count.txt.fillBars"),".png")
knitr::include_graphics(figs_2)
```
## 各样品中亲本差异表达基因数量概述 {#snp-gene-sum}

样品组内亲本间显著差异基因数量(P < $opts{p}, FDR < $opts{q})如下表所示。[数量汇总](snp/gene_sum.txt)  [详细列表](snp/gene.txt)

```{r snp-gene-sum}
table_gene <- read.table("snp/gene_sum.txt", sep="\t", header=T)
knitr::kable(table_gene, caption="各组内亲本差异表达基因数量汇总", booktabs=TRUE)
```


## 差异基因分析 {#snp-gene-de}

!;

my @group;
open DATABASE,"<snp/group_compare.txt";
while (<DATABASE>) {
	chomp;
	push @group,$_;
}
$i=0;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	$i++;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

(ref:snp-$i) 样品内$tmp[0]vs$tmp[1]亲本来源间存在一些表达显著上调或下调的基因 (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示$tmp[0]与$tmp[1]亲本来源相比基因表达的差异分布，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调；(B) 热图展示$tmp[0]与$tmp[1]显著差异基因在每个样品中相对丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示基因延染色体的差异分布情况，对于观察某些成簇表达变化的基因，x轴为基因在染色体上排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，大小为相对丰度，颜色和形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of number of differentially expressed genes, fold change and expression level; (B) Heatmap showing significantlly  differentially expressed genes; (C) Manhattan plot showing differentially expressed genes along the chromosome.
[Volcano plot PDF](snp/vol_gene_$tmp[0]vs$tmp[1].pdf)  [Heatmap PDF](snp/heat_gene_$tmp[0]vs$tmp[1]_sig.pdf)  [Manhattan plot PDF](snp/man_gene_$tmp[0]vs$tmp[1].pdf)

```{r snp-$i, fig.cap="(ref:snp-$i)", out.width="99%"}
figs_2 = paste0("snp/", c("vol_gene_$tmp[0]vs$tmp[1]", "heat_gene_$tmp[0]vs$tmp[1]_sig", "man_gene_$tmp[0]vs$tmp[1]"),".png")
knitr::include_graphics(figs_2)
```
$tmp[0]与$tmp[1]相比不同亲本来源显著差异基因信息。[Enriched TXT](snp/gene_$tmp[0]vs$tmp[1]_enriched.txt)  [Depleted TXT](snp/gene_$tmp[0]vs$tmp[1]_depleted.txt)


(ref:GO-snp-$i) 不同亲本来源$tmp[0]vs$tmp[1]组间存在一些表达显著上调或下调的基因功能富含分析。 (P & FDR < 0.01, 超几何分布)。(A) GO生物学过程分类富集结果；(B) GO分子功能分类富集结果；(C) GO细胞组分分类富集结果。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.01, hypergeometric distribution test). (A) Enriched GO terms involved in biological process; (B) Enriched GO terms involved in molecular function; (C) Enriched GO terms involved in celluar component; (D) Enriched KEEG pathway.
[Enriched GO BP](snp/gene_$tmp[0]vs$tmp[1]_enriched_BP.pdf)  [Depleted GO BP](snp/gene_$tmp[0]vs$tmp[1]_depleted_BP.pdf)
[Enriched GO MF](snp/gene_$tmp[0]vs$tmp[1]_enriched_MF.pdf)  [Depleted GO MF](snp/gene_$tmp[0]vs$tmp[1]_depleted_MF.pdf)
[Enriched GO CC](snp/gene_$tmp[0]vs$tmp[1]_enriched_CC.pdf)  [Depleted GO CC](snp/gene_$tmp[0]vs$tmp[1]_depleted_CC.pdf)
[Enriched KEEG](snp/gene_$tmp[0]vs$tmp[1]_enriched_KK.pdf)  [Depleted KEEG](snp/gene_$tmp[0]vs$tmp[1]_depleted_KK.pdf)
[Enriched GO BP_enrichMap](snp/gene_$tmp[0]vs$tmp[1]_enriched_BP_enrichMap.pdf)  [Depleted GO BP_enrichMap](snp/gene_$tmp[0]vs$tmp[1]_depleted_BP_enrichMap.pdf)
[Enriched GO BP_cnetplot](snp/gene_$tmp[0]vs$tmp[1]_enriched_BP_cnetplot.pdf)  [Depleted GO BP_cnetplot](snp/gene_$tmp[0]vs$tmp[1]_depleted_BP_cnetplot.pdf)
[Enriched GO BP_plotGOgraph](snp/gene_$tmp[0]vs$tmp[1]_enriched_BP_plotGOgraph.pdf)  [Depleted GO BP_plotGOgraph](snp/gene_$tmp[0]vs$tmp[1]_depleted_BP_plotGOgraph.pdf)

```{r GO-snp-$i, fig.cap="(ref:GO-snp-$i)", out.width="49%"}
figs_2 = paste0("snp/gene_$tmp[0]vs$tmp[1]_", c("enriched_BP", "depleted_BP","enriched_MF", "depleted_MF","enriched_CC", "depleted_CC","enriched_KK", "depleted_KK"),".png")
knitr::include_graphics(figs_2)
```

!;
}
}
if ($opts{C} eq "TRUE") {
open OUTPUT,">$opts{o}04-h-cross.Rmd";
print OUTPUT qq!
# 跨物种比较pollen, ovule, stigma间差异基因 {#cross-cross}

实验物种包括拟南芥Arabidopsis thaliana Col和Cvi生态型，以及相近物种Arabidopsis lyrata、Brassica rapa、Capsella rubella；以讨论物种间隔离形成的原因。


## 跨物种比较样品信息 {#cross-sample}

```{r design_cross}
table_design_cross <- read.table("cross/design.txt", sep="\\t", header=T)
knitr::kable(table_design_cross, caption="样品详细及分组信息总结。", booktabs=TRUE)
```

## 样品间或组间Pearson相关分析 {#cross-cross-cor}

(ref:cross-group-cor) 基于各样品组皮尔森相关系数。Pearson correlation of all groups. [PDF](cross/heat_cor_groups.pdf)

```{r cross-group-cor, fig.cap="(ref:cross-group-cor)", out.width="99%"}
knitr::include_graphics("cross/heat_cor_groups.png")
```

(ref:cross-sample-cor) 基于各样品组基因RPM值计算皮尔森相关系数。Pearson correlation of all samples. [PDF](cross/heat_cor_samples.pdf)

```{r cross-sample-cor, fig.cap="(ref:cross-sample-cor)", out.width="99%"}
knitr::include_graphics("cross/heat_cor_samples.png")
```

## 跨物种组织间差异表达基因数量概述 {#cross-gene-sum}

跨物种组织间显著差异基因数量(P < $opts{p}, FDR < $opts{q})如下表所示。[数量汇总](cross/gene_sum.txt)  [详细列表](cross/gene.txt)

```{r cross-gene-sum}
table_gene <- read.table("cross/gene_sum.txt", sep="\t", header=T)
knitr::kable(table_gene, caption="跨物种组织间亲本差异表达基因数量汇总", booktabs=TRUE)
```


## 跨物种组织间差异基因分析 {#cross-gene-de}

!;

my @group;
open DATABASE,"<cross/group_compare.txt";
while (<DATABASE>) {
	chomp;
	push @group,$_;
}
$i=0;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	$i++;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

(ref:cross-$i) 跨物种$tmp[0]vs$tmp[1]间存在一些表达显著上调或下调的基因 (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示$tmp[0]与$tmp[1]相比基因表达的差异分布，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调；(B) 热图展示$tmp[0]与$tmp[1]显著差异基因在每个样品中相对丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示基因沿染色体的差异分布情况，对于观察某些成簇表达变化的基因，x轴为基因在染色体上排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，大小为相对丰度，颜色和形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of number of differentially expressed genes, fold change and expression level; (B) Heatmap showing significantlly  differentially expressed genes; (C) Manhattan plot showing differentially expressed genes along the chromosome.
[Volcano plot PDF](cross/vol_gene_$tmp[0]vs$tmp[1].pdf)  [Heatmap PDF](cross/heat_gene_$tmp[0]vs$tmp[1]_sig.pdf)  [Manhattan plot PDF](cross/man_gene_$tmp[0]vs$tmp[1].pdf)

```{r cross-$i, fig.cap="(ref:cross-$i)", out.width="99%"}
figs_2 = paste0("cross/", c("vol_gene_$tmp[0]vs$tmp[1]", "heat_gene_$tmp[0]vs$tmp[1]_sig", "man_gene_$tmp[0]vs$tmp[1]"),".png")
knitr::include_graphics(figs_2)
```
$tmp[0]与$tmp[1]相比不同物种来源显著差异基因信息。[Enriched TXT](cross/gene_$tmp[0]vs$tmp[1]_enriched.txt)  [Depleted TXT](cross/gene_$tmp[0]vs$tmp[1]_depleted.txt)  [Nosig TXT](cross/gene_$tmp[0]vs$tmp[1]_nosig.txt)  


(ref:GO-cross-$i) 不同物种来源$tmp[0]vs$tmp[1]组间存在一些表达显著上调或下调的基因功能富含分析。 (P & FDR < 0.01, 超几何分布)。(A) GO生物学过程分类富集结果；(B) GO分子功能分类富集结果；(C) GO细胞组分分类富集结果。
Compared with $tmp[1], $tmp[0] are enriched and depleted for certain genes. (P & FDR < 0.01, hypergeometric distribution test). (A) Enriched GO terms involved in biological process; (B) Enriched GO terms involved in molecular function; (C) Enriched GO terms involved in celluar component; (D) Enriched KEEG pathway.
[Enriched GO BP](cross/gene_$tmp[0]vs$tmp[1]_enriched_BP.pdf)  [Depleted GO BP](cross/gene_$tmp[0]vs$tmp[1]_depleted_BP.pdf)
[Enriched GO MF](cross/gene_$tmp[0]vs$tmp[1]_enriched_MF.pdf)  [Depleted GO MF](cross/gene_$tmp[0]vs$tmp[1]_depleted_MF.pdf)
[Enriched GO CC](cross/gene_$tmp[0]vs$tmp[1]_enriched_CC.pdf)  [Depleted GO CC](cross/gene_$tmp[0]vs$tmp[1]_depleted_CC.pdf)
[Enriched KEEG](cross/gene_$tmp[0]vs$tmp[1]_enriched_KK.pdf)  [Depleted KEEG](cross/gene_$tmp[0]vs$tmp[1]_depleted_KK.pdf)
[Enriched GO BP_enrichMap](cross/gene_$tmp[0]vs$tmp[1]_enriched_BP_enrichMap.pdf)  [Depleted GO BP_enrichMap](cross/gene_$tmp[0]vs$tmp[1]_depleted_BP_enrichMap.pdf)
[Enriched GO BP_cnetplot](cross/gene_$tmp[0]vs$tmp[1]_enriched_BP_cnetplot.pdf)  [Depleted GO BP_cnetplot](cross/gene_$tmp[0]vs$tmp[1]_depleted_BP_cnetplot.pdf)
[Enriched GO BP_plotGOgraph](cross/gene_$tmp[0]vs$tmp[1]_enriched_BP_plotGOgraph.pdf)  [Depleted GO BP_plotGOgraph](cross/gene_$tmp[0]vs$tmp[1]_depleted_BP_plotGOgraph.pdf)

```{r GO-cross-$i, fig.cap="(ref:GO-cross-$i)", out.width="49%"}
figs_2 = paste0("cross/gene_$tmp[0]vs$tmp[1]_", c("enriched_BP", "depleted_BP","enriched_MF", "depleted_MF","enriched_CC", "depleted_CC","enriched_KK", "depleted_KK"),".png")
knitr::include_graphics(figs_2)
```

!;

}

print OUTPUT qq!
## 比较组间共有和特有基因 {#cross-venn}

!;

## 读group venn文件
open DATABASE,"<cross/group_venn.txt";
while (<DATABASE>) {
	chomp;
	push @venn,$_;
}
close DATABASE;

my $i=0;
foreach (@venn) {
	chomp;
	$i++;
	my @tmp=split/\t/; # sampleA and sampleB
	my $venn_list2;
	foreach $tmp (@tmp) {
		$venn_list2.=$tmp;
	}
	$tmp[2]="C" unless defined($tmp[2]);
	$tmp[3]="D" unless defined($tmp[3]);
	$tmp[4]="E" unless defined($tmp[4]);
	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";

print OUTPUT qq!

### $venn_list

(ref:cross-venn-$i) 维恩图展示各比较组差异基因间的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [$venn_list venn PDF](cross/gene.txt.venn$venn_list.pdf)  [$venn_list venn XLS](cross/gene.txt.venn$venn_list2.xls)  


```{r cross-venn-$i, fig.cap="(ref:cross-venn-$i)", out.width="99%"}
figs_2 = paste0("cross/gene.txt.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}
close OUTPUT;

}



open OUTPUT,">$opts{o}05-reference.Rmd";
print OUTPUT qq!
# Reference {#result-reference}

!;



###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

if ($opts{e} eq "TRUE") {
	`Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"`;
}

open ACCESS,">$opts{V}/.htaccess";
print ACCESS qq!
# AuthName must have, "" not allow blank, not need change
AuthName "Need user name and password"
AuthType Basic
AuthUserFile /mnt/bai/yongxin/bin/config/users
require valid-user
!;

`chmod +x $opts{V}/.htaccess `;
`rm -f $opts{V}/$opts{V}`;
print "Result please visiting http://bailab.genetics.ac.cn/report/rna/$opts{V}\n";
