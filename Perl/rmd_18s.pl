#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    rmd_16s.pl -i inpute_file -o output_file -d design.txt -l library.txt -c group_compare.txt -v group_venn.txt -g groupID
Function: write 16S report in Rbookdown
Command:  -i inpute file name (Must)
          -o output file directory (Must)
          -e TRUE, executate Rscript
          -d design.txt
          -g group1
          -D group2
          -l library.txt
          -c group_compare_elite.txt
          -v group_venn.txt
          -s summary.txt
          -t group_tern.txt
          -a abundance threshold, default 0.001
          -n name of user
          -b version for output directory
          -S TRUE or FLASE, whether report simplified report, default FALSE
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.3
Update:   2017/6/26
Notes:    1.0 outpu html report and logo
          1.1 add reference bib, replace reference to link bib
          1.2 add summary.txt for project information; add check whether have CPCoA and venn; Alpha change to raw_OTU table and add rarefaction plot
          1.3 add ternery plot, simplify report table and figures
          1.4 add venn number, split into each figures
          1.5 add phylum+class taxonomy; add manhattan plot in order
\n!
    )
}
###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'a:i:o:d:e:l:c:v:n:h:s:g:D:t:b:S:', \%opts );
$opts{a}=0.001 unless defined($opts{a});
$opts{h}=0 unless defined($opts{h});
$opts{o}="./" unless defined($opts{o}); # work directory
$opts{s}="doc/summary.txt" unless defined($opts{s});
$opts{d}="doc/design.txt" unless defined($opts{d});
$opts{l}="doc/library.txt" unless defined($opts{l});
$opts{c}="doc/group_compare_elite.txt" unless defined($opts{c});
$opts{v}="doc/group_venn.txt" unless defined($opts{v});
$opts{t}="doc/group_tern.txt" unless defined($opts{t});
$opts{e}="TRUE" unless defined($opts{e}); # report elite version
#$opts{b}=$opts{n} unless defined($opts{b}); # output directory, default same with username
$opts{g}="genotype" unless defined($opts{g}); # default group column name
$opts{S}="FALSE" unless defined($opts{S}); # report elite version
$opts{c}="$opts{c}.elite" if $opts{S};
$opts{c}=~s/\.elite//g unless -e "$opts{c}";

&usage unless (exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

# Prepare relative files
print "Clean report enviroment and prepare relative files\n";
`rm -f *.Rmd`;
`rm -fr _bookdown_files`;
`rm -fr $opts{b}`;
`cp -f -r /mnt/bai/yongxin/ref/amplicon/rmd/* $opts{o}`;
`sed -i 's/html/$opts{b}/g' _bookdown.yml`;

my %list;
open LIST,"<$opts{s}";
while (<LIST>) {
	chomp;
	my @tmp=split/\t/;
	$list{$tmp[0]}=$tmp[1];
}
close LIST;

open OUTPUT,">$opts{o}index.Rmd";
print OUTPUT qq!--- 
title: "$list{title}"
author:
- 实验数据：$list{client}
- 生信分析：中国科学院遗传与发育生物学研究所白洋组[刘永鑫](http://blog.csdn.net/woodcorpse)
- 联系方式：yxliu\@genetics.ac.cn 010-64808722
- 项目编号：$list{project}
- 项目周期：$list{period}
- 官方网站：http://bailab.genetics.ac.cn/
- 公 众 号：[植物微生物组](http://mp.weixin.qq.com/s/QkgNlzK_rpauKzSs2fUd7A)  [宏基因组](http://mp.weixin.qq.com/s/bPdCViMURAByLBcGqOmE4A)
date: '`r Sys.Date()`'
documentclass: article
bibliography: [16s.bib]
link-citations: yes
biblio-style: apalike
---

```{r setup, include=FALSE}
library(knitr)
output = opts_knit\$get("rmarkdown.pandoc.to")
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
mtime = function(files){
  lapply(Sys.glob(files), function(x) file.info(x)\$mtime)
}
set.seed(718)
```

```{asis, echo=html}
# Bailab, SKLPG/CEPAMS, IGDB, CAS {-}
```

```{r cover, eval=html, out.width="99%"}
#figs_1 = paste0("figure/slide", c("1", "2"),"_raw.jpg")
#knitr::include_graphics(figs_1)
knitr::include_graphics("figure/banner.png")
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
#`sed ':a;N;\$!ba;s/^/ID\\t/g' doc/design.txt > doc/design.xls`;
`ln doc/design.txt doc/design.xls -f`;
print OUTPUT qq!
# 课题设计 {#project_design}

$list{design}样品信息如 Table \\\@ref(tab:design) 所示，[完整实验设计表design.xls](doc/design.xls)。

```{r design}
table_design = read.table("doc/design.txt", sep="\\t", header=T)
table_design=head(table_design,n=10)
knitr::kable(table_design, caption="样品信息前10行。Top 10 lines of sample information.", booktabs=TRUE)
```
!;
close OUTPUT;

if ($opts{S} eq "FALSE") {
open OUTPUT,">$opts{o}03-scheme.Rmd";
print OUTPUT qq!
# 课题方案 {#project_scheme}

我的材料有那些菌？taxonomy tree, phylogenetic tree.  

实验组和对照组间是否存在不同？alpha diversity, beta diversity.  

具体有那些不同？Differentially abundance taxonomy and OTU.  

整个分析流程包含以下10部分内容：报告测序数据质控；测序数据过滤及各步骤统计、样品数据量和长度分布；Alpha多样性分析: Shannon entropy和observed OTU；Beta多样性分析: 采用bray curtis和weighted unifrac距离计算距离的主坐标轴分析(PCoA/MDS)；限制条件的PCoA分析(CPCoA/CCA/RDA); 分类树及进化树展示OTU物种信息及进化关系；各分类级别丰度分析：包括门、纲、目、科、属水平；差异OTU分析：包括火山图、热图、曼哈顿图展示差异OTU数量、丰度、变化样式及分类学信息；组间差异OTU比较，观察不同组间的分类学样式，以及共有或特有OTU；其它有待进一步分析的内容，如OTU调控网络构建等。  

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


(ref:scheme-read-fastqc) 测序Reads质量评估。HiSeq2500产出Clean reads左端(A)和右端(B)各250 bp数据质量评估，选取测序reads碱基质量分布判断建库或测序有无异常。双端数据raw和clean reads左端(C)和右端(D)接头及引物污染情况分布，接头去除干净与否、及有效数据比例评估。  Quality control of raw reads [\@andrews2010fastqc]

```{r scheme-read-fastqc, fig.cap="(ref:scheme-read-fastqc)"}
knitr::include_graphics("figure/fig1.fastqc.png")
```

**2. 样品提取及过滤各步骤统计；Statistics of reds filter processes**

(ref:scheme-read-summary) 统计文库处理过程及样品可用数据量。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B). 柱状图展示各样品的数据量分布，最小值也大于2万，大部分在12万左右，完全符合实验设计要求；(C) 可用数据的长度分布，可以看到本实验扩增子长度范围集中在360-390 bp，主峰位于370-380 bp间。  Statistics of reads filter processes in libraries and data size of samples. (A) Bar plot showing reads count of each library in read filter process; (B) Bar plot showing reads counts of each sample; (C) Length distribution of amplicons [\@caporaso2010qiime, \@edgar2013uparse].

```{r scheme-read-summary, fig.cap="(ref:scheme-read-summary)"}
knitr::include_graphics("figure/fig2.summary.png")
```

**3. Alpha多样性分析；Alpha (α) diversity**

(ref:scheme-sample-alpha) Alpha多样性展示各组间微生物多样性，方法采用(A) Shannon index，包括样品的可操作分类单元(operational taxonomic unit, OTU)数量及种类丰度信息；(B) Observed OTUs index，只包括样品OTU种类信息。图中KO(knock out)代表基因敲除突变体，OE(overexpression)代表基因过表达株系，WT(wild-type)代表野生型。附表有各种间t-test方法统计的p-value水平。此外还可计算chao1和PD whole tree等方法下的多样性分析。[各Alpha多样性计算方法详细](http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html)  Within sample diversity (α-diversity) measurements among each genotype. (A) Shannon index, estimated species richness and evenness; (B) Observed OTUs index, only calculate species richness. These results indicate genotype not significantly change microbial diversity. The horizontal bars within boxes represent median. The tops and bottoms of boxes represent 75th and 25th quartiles, respectively. The upper and lower whiskers extend 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. All outliers are plotted as individual points [\@edwards2015structure].

```{r scheme-sample-alpha, fig.cap="(ref:scheme-sample-alpha)"}
knitr::include_graphics("figure/fig3.alpha.png")
```

**4. Beta多样性分析；Beta (β) diversity **

(ref:scheme-sample-beta) 采用主坐标轴分析展示第1/2坐标轴下各组间微生物组差异(dissimilarity)，距离计算方法采用(A) bray curtis; (B) weighted unifrac. 如图A中可以看到坐标轴1可以解释24.15%的变异，坐标轴2可以解释12.32%的变异，KO与WT较为相似；而OE在第一轴方向上明显与WT分开，表明其微生物组呈现明显变化；同时还发现KO1中存在三个样品存在明显异常。  Principal coordinate analysis (PCoA) using the (A) bray curtis metric and (B) weighted unifrac metric shows dissimilarity of microbial communities. The result indicates that the largest separation is between WT and OE (PCoA 1) and the second largest source of variation is between WT and KO (PCoA 2) [\@edwards2015structure].

```{r scheme-sample-beta, fig.cap="(ref:scheme-sample-beta)"}
knitr::include_graphics("figure/fig4.beta.png")
```

**5. 限制条件下的主坐标轴分析；Constrained principal coordinate analysis**

(ref:scheme-sample-CPCoA) 以基因型为条件分析贡献率及组间差异；分析表明基因型可解释微生物组的22.7%的变异，且各基因型间均可明显分开，且KO和OE的重复又能很好聚在一起，表明不同基因对微生物组的群落结构有明显的调控作用，且不同突变体和过表达株系的位点和生物学重复间表现出良好的可重复性。  Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype (22.7% of the overall variance; p < 0.05) [\@bulgarelli2015structure].

```{r scheme-sample-CPCoA, fig.cap="(ref:scheme-sample-CPCoA)"}
knitr::include_graphics("figure/fig5.CPCoA.png")
```

**6. 分类树及进化树展示OTU物种信息及进化关系；Taxonomy and phylogenetic tree of OTU**

(ref:scheme-sample-tree) 样品中高丰度(>0.5%)OTU的分类树和系统发生学分析。(A)分类树，其中OTU按分类学的科水平进行背景高亮着色，显示本研究中主要丰度的细菌科；(B)系统发生树，按门水平进行着色，结果表明细菌的物种注释信息与16S的序列发生树的进化关系高度一致。  Taxonomy and phylogenetic tress show high abundance OTU (>0.5%), and their family and phylum annotation of taxonomy [\@asnicar2015compact, \@yu2017ggtree]. 

```{r scheme-sample-tree, fig.cap="(ref:scheme-sample-tree)"}
knitr::include_graphics("figure/fig6.tree.png")
```

**7. 分类学不同分类级别的丰度分析；Differentially abundance of bacterial in each taxonomy level**

(ref:scheme-sample-tax) 柱状图展示各类微生物组分类学门水平相对丰度。(A) 堆叠柱状图，X轴为各样品组，Y轴为各门类相对百分比，只列出了丰度大于0.1%的门，其它所有门归入Low Abundance类。(B). 条形图展示最高丰度的五大菌门平均丰度及标准误，我们可以观察到与WT相比，各基因型的Proteobacteria丰度降低，而Actinobacteria丰度升高。注: 分类学注释可从门、纲、目、科、属五个级别进行丰度可视化及差异统计分析。  Bar plot showing phyla abundances in each genotype. (A). Stack plot showing high abundance (>0.1%) phylum; (B). Bar plot showing top 5 phylum abundance and error bar in each genotype. All the KO and OE were show enriched in Actinobacteria and depleted in Proteobacteria. Note: Differentially abundance taxonomy can analyze under phylum, order, class, family and genus level [\@bulgarelli2015structure, \@lebeis2015salicylic].

```{r scheme-sample-tax, fig.cap="(ref:scheme-sample-tax)"}
knitr::include_graphics("figure/fig7.tax.png")
```

**8. 差异OTUs分析；Differentially abundance OTUs**

(ref:scheme-sample-otu) KO1基因型存在一些丰度显著上调或下调的OTU (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示KO与WT相比OTU的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调，图中数字代表显著差异OTU数量，形状代表OTU的门水平物种注释；（B）热图展示KO与WT显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低；可以看到我们找到的差异OTU在每组样品中重复非常好，同时也发现了在beta diversity分析中发现的KO1中存在的两个异常样品应该为KO1.7, KO1.8, 需检查实验材料准备了取材步骤有无问题？或补弃样品重复（C）曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调。  KO1 are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs of KO1 compared WT; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs. These results show Actinobacterial has more enriched OTUs [\@bai2015functional, \@edwards2015structure, \@zgadzaj2016root].

```{r scheme-sample-otu, fig.cap="(ref:scheme-sample-otu)"}
knitr::include_graphics("figure/fig8.otu.png")
```

**9. 组间差异OTU比较；Compare differentially abundance OTUs among groups**

(ref:scheme-sample-overlap) 比较组间差异OTU的分类学样式、共有或特有。(A) 饼形图展示各种差异OTU细菌门水平分类比例。中间数字为所有显著差异OTU的数目。可以看到KO1与KO2样式相似，OE1与OE2样式相似。且上调OTU较多为Actinobacteria，而下调OTU绝大多数为Proteobacteria。(B) 维恩图展示各基因型差异OTUs间的共有和特有数量。图中所显各基因型组间重复间大部分OTUs共有；而且还发现KO和OE还存在一些相似变化样式的OTUs。  Taxonomy, common and unique OTUs in each group. (A) Pie charts show phyla of bacterial OTUs identified as either enriched or depleted in each genotype compared with WT. The number of OTUs in each category is noted inside each donut. (B) Venn diagrams show common and unique OTUs in each group [\@lebeis2015salicylic].

```{r scheme-sample-overlap, fig.cap="(ref:scheme-sample-overlap)"}
knitr::include_graphics("figure/fig9.overlap.png")
```

**10. 其它数据分析过程中发现的有意思的点，商讨后，有意义的深入分析；Other points and ideas for further discussion and analysis **
!;
close OUTPUT;
}



## 读取文库列表文件
open DATABASE,"<$opts{l}";
<DATABASE>;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @library,$tmp[0] if defined($tmp[0]); # 判断是否非空的行
}
close DATABASE;

open OUTPUT,">$opts{o}04-a-sequenceQuality.Rmd";
print OUTPUT qq!
# 测序质量总结 {#sequencing_quality_summary}

## 测序质量评估 {#sub-sequence-qc}

说明：16S扩增子测序数据主要来自HiSeq2500产出的双端各250 bp (PE250)数据，优势是读长长且价格适中(性价比高)。HiSeqX PE150和MiSeq PE300也比较常见，但PE150过短分辨率低，而PE300价格高且末端序列质量过低。此外454在之前研究较多且设备已经停产，PacBio读长长可直接测序16S全长1.5kb代表未来的趋势。测序公司通常会返回raw data和clean data两种数据，raw data为测序获得的原始数据，而clean data则为去除含有接头序列及测序不确定N比例较高的结果，通常直接采用clean data进行质量评估及后续分析。数据质量评估结果中测序reads碱基质量分布图，常用于判断建库或测序有无异常。序列重复情况分布，判断原始序列的DNA质量、重复序列比例及PCR扩增重复情况，如重复序列较高可能某些菌高丰度或PCR扩增导致，对低丰度菌的结果影响较大。

!;
close OUTPUT;

foreach $library (@library) {
open OUTPUT,">>$opts{o}04-a-sequenceQuality.Rmd";
`cp -r clean_data/${library}_*.html result/`;
`cp -r clean_data/${library}_1_fastqc/ result/`;
`cp -r clean_data/${library}_2_fastqc/ result/`;
print OUTPUT qq!
### 文库${library}质量评估

(ref:quality-fastqc-${library}) 测序Reads质量评估文库${library}。Clean reads左端(A)和右端(B)数据质量评估；clean reads左端(C)和右端(D)序列重复情况分布。Quality control of clean reads [HTML report of library ${library}_1](result/${library}_1_fastqc.html)  [HTML report of library ${library}_2](result/${library}_2_fastqc.html)

```{r quality-fastqc-${library}, fig.cap="(ref:quality-fastqc-${library})", out.width="49%"}
figs_1 = paste0("result/${library}_", c("1_fastqc/Images/per_base_quality", "2_fastqc/Images/per_base_quality", "1_fastqc/Images/duplication_levels", "2_fastqc/Images/duplication_levels"),".png")
knitr::include_graphics(figs_1)
```

(ref:quality-split-${library}) 文库${library}各样品按barcode拆分获得的高质量测序数据，按实验设计组着色。Distribution of sequenced reads of samples in library ${library}. Samples were colored by group information. 1 Million = 10^6^. [PDF](result/stat_lib_split_${library}.pdf)

```{r quality-split-${library}, fig.cap="(ref:quality-split-${library})"}
knitr::include_graphics("result/stat_lib_split_${library}.png")
```

!;
close OUTPUT;
}

open OUTPUT,">>$opts{o}04-a-sequenceQuality.Rmd";
`ln result/qc.sum result/qc.xls -f`;
`echo -e "Library\tCount\tLength">temp/temp`;
`cat temp/temp result/length.txt > result/length.xls`;
`rm temp/temp`;
print OUTPUT qq!
## 测序数据预处理总结 {#sub-sequence-summary}

(ref:quality-sum) 测序文库数据量和长度分布。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B) 折线图展示各测序文库中序列的长度分析。Data size and length distribution of sequencing libraries. (A) Bar plot showing reads count of each library in read filter process. (B) Line plot showing reads length distribution of each library. [Sum PDF](result/stat_lib_qc_sum.pdf) [Sum XLS](result/qc.xls) [Length PDF](result/stat_lib_length.pdf) [Length XLS](result/length.xls)

```{r quality-sum, fig.cap="(ref:quality-sum)"}
knitr::include_graphics(c("result/stat_lib_qc_sum.png","result/stat_lib_length.png"))
```

!;
close OUTPUT;



open OUTPUT,">$opts{o}04-b-tree.Rmd";
print OUTPUT qq!
# OTU表统计、分类和进化树分析 {#result-tree}

## OTU表基本信息 {#sub-result-table}

```{r otutable}
tab = read.table("result_k1-c/otu_table.txt", sep="\\t", header=T)
tab=head(tab,n=10)
knitr::kable(tab, caption="OTU表前10行。Top 10 lines of OTU table.", booktabs=TRUE)
```

原始数据经过去冗余、过滤出现1次的序列、按97%聚类、自身比对去嵌合体、与RDP数据比对去嵌合体、与Green gene比对去除非细菌序列。原始含taxonomy注释OTU表下载[BIOM](result/otu_table_tax.biom) [TXT](result/otu_table_tax.txt)。基本统计信息如下：

```{r otutable-sum}
tab = read.table("result/otu_table.sum", sep=":", header=F)
colnames(tab)=c("Types","Values")
tab=head(tab,n=10)
knitr::kable(tab, caption="OTU表基本统计信息。Summary of OTU table.", booktabs=TRUE)
```

原始OTU表只保留在至少在一个样品中存在丰度大于$opts{a}的OTU，并去除受Cyanobacteria和Chloroflexi菌门细菌后，接下来的数据进行差异统计分析。原始OTU下载 [BIOM](result_k1-c/otu_table_tax.biom) [TXT](result_k1-c/otu_table_tax.txt)。基本统计信息如下：

```{r otutable-filter-sum}
tab = read.table("result_k1-c/otu_table.sum", sep=":", header=F)
colnames(tab)=c("Types","Values")
tab=head(tab,n=10)
knitr::kable(tab, caption="OTU表基本统计信息。Summary of OTU table.", booktabs=TRUE)
```

## 高丰度OTU物种分类树分析 {#sub-result-graphlan}

(ref:tree-graphlan) 高丰度OTU(>0.5%)物种注释分类树，按分类学目(A. order)、科(B. family)、属(C. genus)水平进行文本标签注释，结果可以看到本实验中鉴定的细菌OTU主要分布于不同分类级别的哪些目、科、属。Taxonomy tress show high abundance OTU (>0.5%), and their order, family and genus annotaion of taxonomy. [order PDF](result/tax_order.pdf)  [family PDF](result/tax_family.pdf)  [genus PDF](result/tax_genus.pdf)

```{r tree-graphlan, fig.cap="(ref:tree-graphlan)"}
if ($opts{S}) {figs_2 = paste0("result/tax_", c("order"),".png")
}else{figs_2 = paste0("result/tax_", c("order", "family", "genus"),".png")}
knitr::include_graphics(figs_2)
```

## 高丰度OTU进化分析 {#sub-result-ggtree}

(ref:tree-ggtree) 高丰度OTU系统发生树分析(>0.5%)，按分类学门(A. phylum)、纲(B. class)、目(C. order)水平进行着色，结果可以看到本实验中鉴定的细菌OTU主要分布于那些分类级别，同时表明细菌的物种注释信息与16S的序列发生树的进化关系高度一致。Phylogenetic tress show high abundance OTU (>0.5%), and their phylum, class and order annotaion of taxonomy. [phylum PDF](result/ggtree_phylum.pdf); [class PDF](result/ggtree_class.pdf); [order PDF](result/ggtree_order.pdf).

```{r tree-ggtree, fig.cap="(ref:tree-ggtree)", out.width="99%"}
if ($opts{S}) {figs_1 = paste0("result/ggtree_", c("phylum"),".png")
}else{figs_1 = paste0("result/ggtree_", c("phylum", "class", "order"),".png")}
knitr::include_graphics(figs_1)
```
!;
close OUTPUT;


open OUTPUT,">$opts{o}04-d-diversity.Rmd";
`ln result/alpha.txt result/alpha.xls -f`;
`ln result_k1-c/beta.txt result_k1-c/beta.xls -f`;
print OUTPUT qq!
# 样品与组多样性分析 {#result-diversity}

## 样品与组间差异分析 {#result-diversity-cor}

(ref:group-cor) 基于各样品组OTU均值计算差异(1-皮尔森相关系数)。Dissimilarity (1 - Pearson correlation) of all groups. [PDF](result_k1-c/heat_cor_groups.pdf)

```{r group-cor, fig.cap="(ref:group-cor)", out.width="99%"}
knitr::include_graphics("result_k1-c/heat_cor_groups.png")
```

(ref:sample-cor) 基于各样品组OTU相对丰度计算差异(1-皮尔森相关系数)。Dissimilarity (1 - Pearson correlation) of all samples. [PDF](result_k1-c/heat_cor_samples.pdf)

```{r sample-cor, fig.cap="(ref:sample-cor)", out.width="99%"}
knitr::include_graphics("result_k1-c/heat_cor_samples.png")
```

## Alpha多样性分析：物种丰富度、均匀度 {#result-diversity-alpha}

### 各样品Alpha多样性数值 {#result-diversity-alpha-value}

各样品常用四种Alpha多样性计算方法结果见如 Table \\\@ref(tab:alpha) 所示。[完整表格下载TXT](result/alpha.xls)

```{r alpha}
table_alpha = read.table("result/alpha.xls", sep="\\t", header=T, row.names = 1)
table_alpha=head(table_alpha,n=10)
knitr::kable(table_alpha, caption="样品四种Alpha多样性结果前10行", booktabs=TRUE)
```

### 各样品组Alpha多样性分布图 {#result-diversity-alpha-figure}

[QIIME绘制稀释曲线 Rarefraction plots](result/rarefaction_plots.html)

(ref:div-alpha) 箱线图展示各样品及组的微生物组Alpha多样性，方法采用(A) Shannon index，包括样品的可操作分类单元(operational taxonomic unit, OTU)种类(richness)及丰度(evenness)信息；(B) Observed OTUs index，只包括样品OTU种类信息。(C) Chao1 index,基于样品测序中单拷贝OTU(饱合情况)估算样品物种种类的方法; (D) PD whole tree index, 多样性评估时考虑OTU间的进化关系，通常进化关系相近的物种可能存在丰度更相关。图中KO(knock out)代表基因敲除突变体，OE(overexpression)代表基因过表达株系，WT(wild-type)代表野生型。附文本有t-test方法统计各组间是否存在显著差异的p-value水平。
[Shannon TXT](result/alpha_shannon_stats.txt)  [observed_otus TXT](result/alpha_observed_otus_stats.txt)  [chao1 TXT](result/alpha_chao1_stats.txt)  [PD_whole_tree TXT](result/alpha_PD_whole_tree_stats.txt)
Within sample diversity (α-diversity) measurements among each genotype. (A) Shannon index, estimated species richness and evenness; (B) Observed OTUs index, only calculate species richness; (C) Chao1 index, calculate richness based on observed, singletons and doubletons; (D) PD whole tree index, diversity considered the evolution distance as weighted. These results indicate genotype not significantly change microbial diversity. The horizontal bars within boxes represent median. The tops and bottoms of boxes represent 75th and 25th quartiles, respectively. The upper and lower whiskers extend 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. All outliers are plotted as individual points (Edwards et al., 2015).
 [Shannon PDF](result/alpha_shannon.pdf)  [observed_otus PDF](result/alpha_observed_otus.pdf)  [chao1 PDF](result/alpha_chao1.pdf)  [PD_whole_tree PDF](result/alpha_PD_whole_tree.pdf) 

```{r div-alpha, fig.cap="(ref:div-alpha)", out.width="99%"}
#figs_2 = paste0("result/alpha_", c("shannon", "observed_otus", "chao1", "PD_whole_tree"),".png")
if ($opts{S}) {figs_2 = paste0("result/alpha_", c("shannon", "observed_otus"),".png")
}else{figs_2 = paste0("result/alpha_", c("shannon", "observed_otus", "chao1", "PD_whole_tree"),".png")}
knitr::include_graphics(figs_2)
```


## Beta多样性分析：基于排序及降维的主坐标轴分析展示样品及组间的差异

(ref:div-beta) 主坐标轴分析(PCoA)展示第1/2坐标轴下各样品间微生物组差异(dissimilarity)，距离计算方法采用(A) bray curtis; (B) unweighted unifrac; (C) weighted unifrac。[采用Adonis统计各样品组间的显著性差异P值](result_k1-c/beta.xls)。
Principal coordinate analysis (PCoA) using the (A) bray curtis metric, (B) unweighted unifrac metric and (C) weighted unifrac metric shows dissimilarity of microbial communities. [bray_curtis PDF](result_k1-c/beta_pcoa_bray_curtis.pdf)  [unweighted_unifrac PDF](result_k1-c/beta_pcoa_unweighted_unifrac.pdf)  [weighted_unifrac PDF](result_k1-c/beta_pcoa_weighted_unifrac.pdf)  [bray_curtis with sampleID lables PDF](result_k1-c/beta_pcoa_bray_curtis_lab.pdf)

```{r div-beta, fig.cap="(ref:div-beta)", out.width="99%"}
#figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis", "unweighted_unifrac", "weighted_unifrac"),".png")
if ($opts{S}) {figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis"),".png")
}else{figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis", "unweighted_unifrac", "weighted_unifrac"),".png")}
knitr::include_graphics(figs_2)
```
!;

# 检测限制条件主成分分析是否有结果图，没有则不输出
$file = "result_k1-c/CPCoA_$opts{g}.png";
if (-e $file) {
print OUTPUT qq!
## 限制条件下的主坐标轴分析(全部OTU)

(ref:div-CPCoA-a) 以基因型为条件分析其贡献率和样品组间差异。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result/CPCoA_$opts{g}.pdf)  [PDF labels](result/CPCoA_$opts{g}_lab.pdf)  

```{r div-CPCoA-a, fig.cap="(ref:div-CPCoA-a)", out.width="99%"}
knitr::include_graphics("result/CPCoA_$opts{g}.png")
```
	
## 限制条件下的主坐标轴分析(高丰度OTU)

(ref:div-CPCoA) 以基因型为条件分析其贡献率和样品组间差异(筛选至少在一个样品中OTU丰度 > $opts{a} )。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result_k1-c/CPCoA_$opts{g}.pdf)  [PDF labels](result_k1-c/CPCoA_$opts{g}_lab.pdf)  


```{r div-CPCoA, fig.cap="(ref:div-CPCoA)", out.width="99%"}
knitr::include_graphics("result_k1-c/CPCoA_$opts{g}.png")
```

!;


# 检测第二分组是否有多个条件，有则逐个条件输出
$opts{D}=~s/"//g;
@g2_list=split(/,/,$opts{D});
if (@g2_list>1) {
foreach $i(@g2_list) {

$file = "result_k1-c/CPCoA_$opts{g}_${i}.png";
if (-e $file) {
print OUTPUT qq!
### 限制条件下的主坐标轴分析${i}组

(ref:div-CPCoA-${i}) 在${i}组内，以基因型为条件分析其贡献率和样品组间差异。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result_k1-c/CPCoA_$opts{g}_${i}.pdf)  

```{r div-CPCoA-${i}, fig.cap="(ref:div-CPCoA-${i})", out.width="99%"}
knitr::include_graphics("result_k1-c/CPCoA_$opts{g}_${i}.png")
```

!;
}
}
}

}
close OUTPUT;



## 样品组比较
open DATABASE,"<$opts{c}";
while (<DATABASE>) {
	chomp;
	push @group,$_;
}
close DATABASE;
@tax_en=qw#phylum class order family genus#;
@tax_cn=qw#门 纲 目 科 属#;

open OUTPUT,">$opts{o}04-e-taxonomy.Rmd";
print OUTPUT qq!
#	分类级别的差异分析 {#result-taxonomy}

## 各分类级下各比较组间差异数量

样品组在不同分类级别上显著差异taxonomy数量(Pvalue < 0.05, FDR < 0.05)如 Table \\\@ref(tab:taxonomy-sum) 所示。[TXT](result_k1-c/tax_sum.txt)

```{r taxonomy-sum}
table_taxonomy = read.table("result_k1-c/tax_sum.txt", sep="\t", header=T)
knitr::kable(table_taxonomy, caption="样品组显著差异taxonomy数量", booktabs=TRUE)
```

!;

foreach $i (0..4) {

# 开启精简模式，只输入门、目水平差异，纲、科和属跳过
if ($opts{S} eq "TRUE") {
	next if $i==1;
	next if $i==2;
	next if $i==4;
}

#print $tax_cn[$i],$tax_en[$i],"\n";
print OUTPUT qq!
## $tax_cn[$i]水平各组差异分析 {#result-taxonomy-$tax_en[$i]}

(ref:taxonomy-$tax_en[$i]) 柱状图展示各样品组微生物组分类学$tax_cn[$i]水平相对丰度。(A) 堆叠柱状图展示各组平均相对丰度，X轴为各样品组，Y轴为各$tax_cn[$i]类相对百分比，只列出了丰度大于0.1%的$tax_cn[$i]，其它所有$tax_cn[$i]归入Low Abundance类。(B) 条形图展示最高丰度的五大菌$tax_cn[$i]平均丰度及标准误，我们可以观察各样品组$tax_cn[$i]水平上相关丰度的差异及组内生物学重复间的波动范围。
Bar plot showing $tax_en[$i] abundances in each genotype. (A) Stack plot showing high abundance (>0.1%) $tax_en[$i]; (B) Bar plot showing top 5 $tax_en[$i] abundance and error bar in each genotype. [stack PDF](result_k1-c/tax_stack_$tax_en[$i]_top9.pdf)  [bar PDF](result_k1-c/tax_bar_$tax_en[$i]_top5.pdf) [raw Data](result_k1-c/database_$tax_en[$i].pdf)

```{r taxonomy-$tax_en[$i], fig.cap="(ref:taxonomy-$tax_en[$i])", out.width="99%"}
#figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_top9", "bar_$tax_en[$i]_top5"),".png")
if ($opts{S}) {figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_top9"),".png")
}else{figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_top9", "bar_$tax_en[$i]_top5"),".png")}
knitr::include_graphics(figs_2)
```
!;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
$file = "result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.pdf";
#print -e $file,"\n"; # 不存在为uninitialized，存在返回1


if (-e $file) {
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.txt > result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls`;
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.txt > result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls`;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

$tmp[0]与$tmp[1]相比显著差异的分类单元信息如 Table \\\@ref(tab:taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]) 所示。[Enriched TXT](result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls)

```{r taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]}
table_taxonomy_e = read.table("result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
table_taxonomy_d = read.table("result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
table_taxonomy_merge = rbind(table_taxonomy_e,table_taxonomy_d)
#table_taxonomy = table_taxonomy_merge[,1:5]
#table_taxonomy = table_taxonomy_merge
table_taxonomy = head(table_taxonomy_merge,n=10)
knitr::kable(table_taxonomy, caption="样品组$tmp[0] vs $tmp[1]显著差异$tax_cn[$i]前10行；Significantlly different $tax_en[$i].", booktabs=TRUE)
```

(ref:taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1]) 热图展示$tmp[0]vs$tmp[1]在$tax_cn[$i]水平差异分类单元。Heatmap show differentially abundance $tax_en[$i].[PDF](result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.pdf)

```{r taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1], fig.cap="(ref:taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1])", out.width="99%"}
knitr::include_graphics("result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.png")
```

!;
}else{
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.

!;
}
}

}



# phylum + proteobacteria class
print OUTPUT qq!
## 门及变形菌纲水平各组差异分析 {#result-taxonomy-phylumpro}

(ref:taxonomy-phylumpro) 柱状图展示各样品组微生物组分类学门及变形菌纲水平相对丰度。(A) 堆叠柱状图展示各组平均相对丰度，X轴为各样品组，Y轴为各门及变形菌纲类相对百分比，只列出了丰度大于0.1%的门及变形菌纲，其它所有门及变形菌纲归入Low Abundance类。(B) 叠柱状图展示各样品相对丰度，我们可以观察各组内样品间的波动情况，又可以看到组间差异。
Bar plot showing phylumpro abundances in each genotype. (A) Stack plot showing high abundance (>0.1%) phylumpro; (B) Stack plot showing high abundance (>0.1%) phylumpro in each sample. [stack PDF](result_k1-c/tax_stack_phylumpro_abc.pdf)  [stack sample PDF](result_k1-c/tax_stack_phylumpro_sample.pdf)  [raw Data](result_k1-c/database_phylumpro.txt)

```{r taxonomy-phylumpro, fig.cap="(ref:taxonomy-phylumpro)", out.width="99%"}
#figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc", "stack_phylumpro_sample"),".png")
if ($opts{S}) {figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc"),".png")
}else{figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc", "stack_phylumpro_sample"),".png")}
knitr::include_graphics(figs_2)
```
!;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
$file = "result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.pdf";
#print -e $file,"\n"; # 不存在为uninitialized，存在返回1


if (-e $file) {
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.txt > result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls`;
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.txt > result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls`;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

$tmp[0]与$tmp[1]相比显著差异的分类单元信息如 Table \\\@ref(tab:taxonomy-$tmp[0]vs$tmp[1]-phylumpro) 所示。[Enriched TXT](result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls)

```{r taxonomy-$tmp[0]vs$tmp[1]-phylumpro}
table_taxonomy_e = read.table("result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
table_taxonomy_d = read.table("result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
table_taxonomy_merge = rbind(table_taxonomy_e,table_taxonomy_d)
#table_taxonomy = table_taxonomy_merge[,1:5]
#table_taxonomy = table_taxonomy_merge
table_taxonomy = head(table_taxonomy_merge,n=10)
knitr::kable(table_taxonomy, caption="样品组$tmp[0] vs $tmp[1]显著差异门及变形菌纲前10行；Significantlly different phylumpro.", booktabs=TRUE)
```

(ref:taxonomy-phylumpro-$tmp[0]vs$tmp[1]) 热图展示$tmp[0]vs$tmp[1]在门及变形菌纲水平差异分类单元。Heatmap show differentially abundance phylumpro.[PDF](result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.pdf)

```{r taxonomy-phylumpro-$tmp[0]vs$tmp[1], fig.cap="(ref:taxonomy-phylumpro-$tmp[0]vs$tmp[1])", out.width="99%"}
knitr::include_graphics("result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.png")
```

!;
}else{
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.

!;
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
## 比较组间共有和特有目 {#result-order-venn}

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

print OUTPUT qq!

### $venn_list

(ref:order-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/order.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/order.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/order.txt.venn$venn_list2.xls.xls)  


```{r order-venn-$i, fig.cap="(ref:order-venn-$i)", out.width="99%"}
figs_2 = paste0("result_k1-c/order.txt.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}

print OUTPUT qq!
## 比较组间共有和特有科 {#result-family-venn}

!;
#(ref:family-venn) 维恩图展示各基因型差异科间的共有和特有数量。图中所显各基因型组间重复间大部分科共有。 Venn diagrams show common and unique familys in each group.!;
#
#foreach (@venn) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#	$tmp[2]="C" unless defined($tmp[2]);
#	$tmp[3]="D" unless defined($tmp[3]);
#	$tmp[4]="E" unless defined($tmp[4]);
#print OUTPUT qq![$tmp[0]$tmp[1] venn PDF](result_k1-c/family.txt.venn$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4].pdf)  !;
#$venn_list.="\"$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]\"\, ";
#}
#$venn_list=~s/\,\ $//;
#print OUTPUT qq!
#
#```{r family-venn, fig.cap="(ref:family-venn)", out.width="49%"}
#figs_2 = paste0("result_k1-c/family.txt.venn", c($venn_list),".png")
#knitr::include_graphics(figs_2)
#```
#
#!;

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

print OUTPUT qq!

### $venn_list

(ref:family-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/family.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/family.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/family.txt.venn$venn_list2.xls.xls)  


```{r family-venn-$i, fig.cap="(ref:family-venn-$i)", out.width="99%"}
figs_2 = paste0("result_k1-c/family.txt.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}

}

print OUTPUT qq!
## 比较组间差异科分类学样式 {#result-family-pie}

(ref:family-pie) 比较组间差异科的门水平分类学样式。饼形图展示各种差异科细菌门水平分类比例。中间数字为所有显著差异科的数目，第一列为显著上调的科，第二列为显著下调的科，从上到下为各比较组。Pie charts show phylum of bacterial familys identified as either enriched or depleted in each genotype compared with WT. The number of familys in each category is noted inside each donut. !;

foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
print OUTPUT qq!
[$tmp[0]vs$tmp[1] enriched pie PDF](result_k1-c/pie_family_$tmp[0]vs$tmp[1]_enriched.pdf) 
[$tmp[0]vs$tmp[1] depleted pie PDF](result_k1-c/pie_family_$tmp[0]vs$tmp[1]_depleted.pdf) !;
$pie_list.="\"$tmp[0]vs$tmp[1]_enriched\"\, \"$tmp[0]vs$tmp[1]_depleted\"\, ";
}
$pie_list=~s/\,\ $//;

print OUTPUT qq!

```{r family-pie, fig.cap="(ref:family-pie)", out.width="49%"}
figs_2 = paste0("result_k1-c/pie_family_", c(${pie_list}),".png")
knitr::include_graphics(figs_2)
```

!;


close OUTPUT;



open OUTPUT,">$opts{o}04-f-otu.Rmd";
print OUTPUT qq!
#	差异OTUs分析 {#result-otu}

## 各总差异OTUs概述 {#result-otu-sum}

样品组间显著差异OTUs数量(P < 0.05, FDR < 0.05)如 Table \\\@ref(tab:otu-sum) 所示。[TXT](result_k1-c/otu_sum.txt)；

```{r otu-sum}
table_otu = read.table("result_k1-c/otu_sum.txt", sep="\t", header=T)
knitr::kable(table_otu, caption="各样品组间差异OTUs数量汇总", booktabs=TRUE)
```


样品组间显著差异显著差异OTU详细列表(P < 0.05, FDR < 0.05)如 Table \\\@ref(tab:otu) 所示。[TXT](result_k1-c/otu.txt)

```{r otu}
table_otu = read.table("result_k1-c/otu.txt", sep="\t", header=F)
colnames(table_otu) = c("OTU","Sample A vs B","P-value")
table_otu=head(table_otu,n=10)
knitr::kable(table_otu, caption="样品组间显著差异前10个OTU, 完整表下载见上方TXT链接", booktabs=TRUE)
```

[下载OTU各组各样品相对丰度和物种信息 raw Data](result_k1-c/database.txt)

## OTU水平组间差异分析 {#result-otu-da}

!;

foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.txt > result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls`;
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.txt > result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls`;
print OUTPUT qq!
### $tmp[0] vs $tmp[1]

(ref:otu-$tmp[0]vs$tmp[1]) $tmp[0]vs$tmp[1]基因型存在一些丰度显著上调或下调的OTU (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示$tmp[0]与$tmp[1]相比OTU的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调；(B) 热图展示$tmp[0]与$tmp[1]显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调；(D) 曼哈顿图按目水平上色。
$tmp[0] are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs, colored by phylum; (D) Manhattan plot colored by order.
[Volcano plot PDF](result_k1-c/vol_otu_$tmp[0]vs$tmp[1].pdf)  [Heatmap PDF](result_k1-c/heat_otu_$tmp[0]vs$tmp[1]_sig.pdf)  [Manhattan plot phlyum PDF](result_k1-c/man_otu_$tmp[0]vs$tmp[1].pdf) [Manhattan plot order PDF](result_k1-c/man_order_$tmp[0]vs$tmp[1].pdf)

```{r otu-$tmp[0]vs$tmp[1], fig.cap="(ref:otu-$tmp[0]vs$tmp[1])", out.width="99%"}
figs_2 = paste0("result_k1-c/", c("vol_otu_$tmp[0]vs$tmp[1]", "heat_otu_$tmp[0]vs$tmp[1]_sig", "man_otu_$tmp[0]vs$tmp[1]", "man_order_$tmp[0]vs$tmp[1]"),".png")
knitr::include_graphics(figs_2)
```

$tmp[0]与$tmp[1]相比显著差异的OTU如 Table \\\@ref(tab:tab-otu-$tmp[0]vs$tmp[1]) 所示。[Enriched TXT](result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls)

```{r tab-otu-$tmp[0]vs$tmp[1]}
table_otu_e = read.table("result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
table_otu_d = read.table("result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
table_otu_merge = rbind(table_otu_e,table_otu_d)
#table_otu = table_otu_merge[,1:5]
#table_otu = table_otu_merge[order(table_otu_merge[,2]),]
table_otu = head(table_otu_merge,n=10)
knitr::kable(table_otu, caption="样品组$tmp[0] vs $tmp[1]相比显著差异的OTU详细信息；Significantlly different OTU.", booktabs=TRUE)
```

!;
}

print OUTPUT qq!
## 比较组间差异OTU分类学样式 {#result-otu-pie}

(ref:otu-pie) 比较组间差异OTU的分类学样式。饼形图展示各种差异OTU细菌门水平分类比例，图例同上图曼哈顿图。中间数字为所有显著差异OTU的数目，第一列为显著上调的OTU，第二列为显著下调的OTU，从上到下为各比较组。Pie charts show phylum of bacterial OTUs identified as either enriched or depleted in each genotype compared with WT. The number of OTUs in each category is noted inside each donut. !;
my $pie_list;
my $venn_list;
foreach (@group) {
	chomp;
	my @tmp=split/\t/; # sampleA and sampleB
print OUTPUT qq!
[$tmp[0]vs$tmp[1] enriched pie PDF](result_k1-c/pie_otu_$tmp[0]vs$tmp[1]_enriched.pdf) 
[$tmp[0]vs$tmp[1] depleted pie PDF](result_k1-c/pie_otu_$tmp[0]vs$tmp[1]_depleted.pdf) !;
$pie_list.="\"$tmp[0]vs$tmp[1]_enriched\"\, \"$tmp[0]vs$tmp[1]_depleted\"\, ";
}
$pie_list=~s/\,\ $//;

print OUTPUT qq!

```{r otu-pie, fig.cap="(ref:otu-pie)", out.width="49%"}
figs_2 = paste0("result_k1-c/pie_otu_", c(${pie_list}),".png")
knitr::include_graphics(figs_2)
```

!;

if (@venn>=1) {
print OUTPUT qq!
## 比较组间共有和特有OTU {#result-otu-venn}

!;
#(ref:otu-venn) 维恩图展示各基因型差异OTUs间的共有和特有数量。图中所显各基因型组间重复间大部分OTUs共有。 Venn diagrams show common and unique OTUs in #each group.!;
#foreach (@venn) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#	$tmp[2]="C" unless defined($tmp[2]);
#	$tmp[3]="D" unless defined($tmp[3]);
#	$tmp[4]="E" unless defined($tmp[4]);
#print OUTPUT qq![$tmp[0]$tmp[1] venn PDF](result_k1-c/otu.txt.venn$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4].pdf)  !;
#$venn_list.="\"$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]\"\, ";
#}
#$venn_list=~s/\,\ $//;
#print OUTPUT qq!
#
#```{r otu-venn, fig.cap="(ref:otu-venn)", out.width="49%"}
#figs_2 = paste0("result_k1-c/otu.txt.venn", c($venn_list),".png")
#knitr::include_graphics(figs_2)
#```
#
#!;

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

(ref:otu-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/otu.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/otu.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/otu.txt.venn$venn_list2.xls.xls)  


```{r otu-venn-$i, fig.cap="(ref:otu-venn-$i)", out.width="99%"}
figs_2 = paste0("result_k1-c/otu.txt.venn", "$venn_list", ".png")
knitr::include_graphics(figs_2)
```

!;
}

}

#$file = "result_k1-c/ter_sum.txt";
if (-e "result_k1-c/ter_sum.txt") {
print OUTPUT qq!
## 三元图展示OTU

三元图中各部分高亮OTU数量统计信息如 Table \\\@ref(tab:ter-sum) 所示。[Summary table](result_k1-c/ter_sum.txt)

```{r ter-sum}
tab = read.table("result_k1-c/ter_sum.txt", sep="\\t", header=T)
knitr::kable(tab, caption="三元图中各部分高亮OTU数量统计", booktabs=TRUE)
```
!;
open DATABASE,"<$opts{t}";
while (<DATABASE>) {
	chomp;
	@tmp=split/\t/;;
	foreach (@tmp) {
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.txt > result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.xls` if -e "result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.txt";
	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.txt > result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.xls` if -e "result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.txt";
	}
print OUTPUT qq!
### $tmp[0] vs $tmp[1] vs $tmp[2]

(ref:ter-$tmp[0]$tmp[1]$tmp[2]) 三元图展示$tmp[0]、$tmp[1]、和$tmp[2]三组OTU(>1‰)的相对丰度。(A) 各组相对于另外两组显著差异的OTU  [PDF](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2].pdf)。每个点代表一个OTU，位置代表它在各组间的相对比例，大小代表三组的平均丰度。彩色的圆圈代表每组中相对于其它两组显著富集的OTU，位于三角形[左$tmp[0]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[0].xls)、[右$tmp[1]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[1].xls)和[顶$tmp[2]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[2].xls)组所特异显著富集的OTU，分别用绿、橙和红高亮显示；
(B) 三角形左$tmp[0]、右$tmp[1]相对于顶$tmp[2] 组富集的OTU  [PDF](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]venn.pdf)。图中底部[两组共有显著富集OTU用红色显示](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[2]venn.xls)，[左$tmp[0]特异的OTU分别用绿色显示](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[0]venn.xls)，[右组特异显示为橙色](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[1].xls)。
Ternary plot depicting group relative abundances (RAs) of all OTUs (>1‰) for $tmp[0], $tmp[1], and $tmp[2]. (A) Compared with each others, three groups specific OTU. Each point corresponds to an OTU. Its position represents its RA with respect to each group, and its size represents the average across all three groups. Colored circles represent OTUs enriched in one group compared with the others (green in $tmp[0], orange in $tmp[1], and red in $tmp[2] samples), whereas gray circles represent OTUs that are not significantly enriched in a specific groups. (B) Group $tmp[2] (top) as control, $tmp[1] and $tmp[2] specific and common OTU; Colored circles represent OTUs enriched in $tmp[0] (green), $tmp[1] (orange), and common in $tmp[0] and $tmp[1] (red).

```{r ter-$tmp[0]$tmp[1]$tmp[2], fig.cap="(ref:ter-$tmp[0]$tmp[1]$tmp[2])", out.width="99%"}
figs_1 = paste0("result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]", c("", "venn"),".png")
knitr::include_graphics(figs_1)
```

!;
}
}
close OUTPUT;



open OUTPUT,">$opts{o}05-references.Rmd";
print OUTPUT qq!
`r if (knitr:::is_html_output()) '# References {-}'`
!;
close OUTPUT;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

if ($opts{e} eq "TRUE") {
	`Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"`;
}
open ACCESS,">$opts{b}/.htaccess";
print ACCESS qq!
# AuthName must have, "" not allow blank, not need change
AuthName "Need user name and password"
AuthType Basic
AuthUserFile /mnt/bai/yongxin/bin/config/users
require valid-user
!;

`chmod +x $opts{b}/.htaccess `;
`cp figure/banner.png $opts{b}/figure/banner.png`;
`cp -r result/average_plots $opts{b}/result`;

print "Result please visiting http://bailab.genetics.ac.cn/report/16s/$opts{b}\n";
