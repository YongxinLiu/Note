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
Usage:    rmd_16s_train.pl
Function: write 16S train in Rbookdown
Command:  -d design.txt
          -g group1
          -D group2_list, for samples in different batch or condition
          -F group3_list, for samples in different batch or condition
          -l library.txt
          -c group_compare_elite.txt
          -v group_venn.txt
          -s summary.txt
          -t group_tern.txt
          -a abundance threshold, default 0.005
          -n name of user
          -b version for output directory
          -S TRUE or FLASE, whether report simplified report, default FALSE
          -h header line number, default 0
Author:   Liu Yong-Xin, yxliu\@genetics.ac.cn, QQ:42789409
Version:  v1.7
Update:   2017/12/3
Notes:    1.0 outpu html report and logo
\n!
    )
}

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'a:i:o:d:e:l:c:v:n:h:s:g:D:F:t:b:S:', \%opts );
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
$opts{b}="train" unless defined($opts{b}); # output directory, default same with username
$opts{g}="genotype" unless defined($opts{g}); # default group column name
$opts{S}="FALSE" unless defined($opts{S}); # report elite version

&usage unless (exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

my %list;
#open LIST,"<$opts{s}";
#while (<LIST>) {
#	chomp;
#	my @tmp=split/\t/;
#	$list{$tmp[0]}=$tmp[1];
#}
#close LIST;

# Set project default parameter
$list{title}="扩增子16S/ITS分析教程" unless defined($list{title});
$list{partner}="中科院遗传发育所白洋组" unless defined($list{partner});
$list{analyst}="刘永鑫 工程师" unless defined($list{analyst});
$list{contact}="yxliu\@genetics.ac.cn 010-64808722" unless defined($list{contact});
$list{project}="171201" unless defined($list{project});
$list{period}="2017-12-01~ 2018-02-28" unless defined($list{period});
$list{website}="http://bailab.genetics.ac.cn/" unless defined($list{website});
$list{wechat}="[植物微生物组](http://mp.weixin.qq.com/s/07Z0CpLhb6-10Nk7tnej7Q)" unless defined($list{wechat});
$list{logo}="Bailab, SKLPG/CEPAMS, IGDB, CAS {-}" unless defined($list{logo});
$list{logo_file}="bailab" unless defined($list{logo_file});

# Prepare relative files
print "Clean report enviroment and prepare relative files\n";
`rm -f *.Rmd`;
`rm -fr _bookdown_files`;
`rm -fr $opts{b}`;
`cp -f -r /mnt/bai/yongxin/ref/amplicon/rmd/* $opts{o}`;
`mv figure/banner_$list{logo_file}.png figure/banner.png`;
`sed -i 's/html/$opts{b}/g' _bookdown.yml`;


open OUTPUT,">$opts{o}index.Rmd";
print OUTPUT qq!
--- 
title: "$list{title}"
author:
- 制作单位：$list{partner}
- 主 讲 人：$list{analyst}
- 联系方式：$list{contact}
- 项目启动：$list{project}
- 项目周期：$list{period}
- 官方网站：$list{website}
- 公 众 号：$list{wechat}
- 更新日期：
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
# $list{logo}
```

```{r cover, eval=html, out.width="99%"}
knitr::include_graphics("figure/banner.png")
```
!;
close OUTPUT;



open OUTPUT,">$opts{o}01-prepare.Rmd";
print OUTPUT qq!
# 分析前准备 {#prepare}

## 登陆服务器 Login bailab

大家最常用的电脑操作系统有Windows10 x64和Mac OS 10两类，我们提供了这两类用户详细操作指南，请根据自己的系统阅读相应部分即可。如果你的系统不属于以上这两种，由于我们没有进行过测试且不可能花时间测试，不确保本教程可以正常使用。建议使用以上两种系统之一，且最好保持系统更新至最新版。

世界上99%的服务器为Linux操作系统，我们做宏基因组学数据分析也不例外。常用的发行版有Ubuntu 16.04 LTS, CentOS 7。

(ref:pre-ubuntu) Ubuntu桌面登陆效果，但在终端的黑屏下功能更强大

```{r pre-ubuntu, fig.cap="(ref:pre-ubuntu)"}
knitr::include_graphics("http://bailab.genetics.ac.cn/train/fig/1.01_ubuntu.jpg")
```

我们布署了Ubuntu 16.04.3 LTS，特点是系统开源且有基金会支持、半年升级一版、软件安装最方便、甚至有中国定制版。感觉是不是和苹果很像呢？因为他们内核是一样的。

### Windows10用户登陆服务器

Windows用户想要登陆Linux服务器，需要安装支持ssh协议的软件，常用的有Xshell、putty、SSH secure shell等。

(ref:pre-xshell) Xshell终端界面，可多窗口同屏

```{r pre-xshell, fig.cap="(ref:pre-xshell)"}
knitr::include_graphics("http://bailab.genetics.ac.cn/train/fig/1.02_xshell.jpg")
```

这里我们推荐使用功能强大的Xshell(没有Mac版)，黑屏代码飞舞，骇客帝国的即视感！

**1. 下载Xshell5**

- 可以访问[官网www.netsarang.com](http://www.netsarang.com/download/)下载最新版，但需要注册申请家庭/学校免费版，且等待下载链接发送到邮箱；

- [不追求最新版，点我下载之前注册的Xshell-5.0.1337p稳定版本](http://bailab.genetics.ac.cn/share/Xshell-5.0.1337p.exe)。

**2. 安装**

- 双击执行安装，按正常顺序下一步，确定即可(如出现输入注册码，选择Free Home/School即可)；

- 安装完成后，桌面和开始菜单都会出现Xshell的图标。 

**3. 登陆**

- 双击Xshell桌面图标运行程序；

- 点击菜单 `文件` —— `新建`，出现`新建会话属性`窗口；

- 填写服务器基本信息：名称为`bailab`，主机为`210.75.224.110`，确定(如出现会话窗口，选择自己新建的`bailab`，点`连接`即可)；

- 提示输入用户名，请输入自己的用户名——姓名后两个字全拼，如刘永鑫为`yongxin`，王鑫为`wangxin`，并勾选下方的`记住用户名`，点击`确定`；

- (第一次登陆会弹出是否信任服务器的安全提醒，选择是)，等待提示身份验证，输入密码`micro.....`，并勾选下方的`记住密码`，点击`确定`；

如果Xshell安装存在问题或无法正常使用，可以下载另一个登陆Linux软件Putty，推荐下载[绿色版PuTTY 0.70](http://bailab.genetics.ac.cn/train/putty.exe)、[安装版](http://bailab.genetics.ac.cn/train/putty-0.70-installer.msi)或访问[主页下载最新版](https://winscp.net/eng/download.php#putty)。


### Mac 10 用户登陆服务器

Mac本来就是Linux的一个商业定制版，它们的内核相同，所以共享大部分程序和命令。如Linux服务器中的终端程序，在Mac上可是标配。

- 自己查找终端 或 Terminal程序并打开

- 输入如下代码连接服务器，并按提示输入密码，注意在Linux里面输入密码是没有任何显示的，不要以为没打上。

注意：下面的`user`请替换为你自己的用户名，为姓名后两个字全拼，如刘永鑫为`yongxin`，王鑫为`wangxin`
```
ssh user\@210.75.224.110
# (第一次登陆会弹出是否信任服务器的安全提醒，请输入`yes`并回车即可)
# 等待提示输入密码即可


# (选学)Linux小技术，在Mac可能有效

# 每次都输入用户名和IP地址太麻烦，可以简化吗
vi ~/.ssh/config # 按i进入编辑模式，添加如下内容至此文件，注意小写的user为自己的用户名

Host bailab
    HostName 210.75.110
    User user
    Port 22

# 按ESC，再按：切换命令模式，输入wq回车实现保存文件并退出

ssh bailab # 相当于打上面的用户名和IP地址，简单吧，再输入密码就登陆成功了


# 每次都要输密码，密码又长又看不见，可以使用我的电脑免密安全登陆吗？当然可以
ssh-keygen # 产生本地公钥和私钥，有提示按说明按y确认或回车即可

# 下用有2个方法，任选其一：
# 1. 专用命令添加
ssh-copy-id bailab
# 2. 用cat id_rsa.pub显示，然后复制到远程服务器的~/.ssh/authorized_keys中
```
参考：http://mp.weixin.qq.com/s/wTM8J9zVEdl1PcGpJQbeEw

## 欢迎进入Linux世界

现在你已经成功登陆到服务器了，并显示如下欢迎信息：

```
Connecting to 210.75.224.110:22...
Connection established.
To escape to local shell, press 'Ctrl+Alt+]'.

Welcome to Ubuntu 16.04.3 LTS (GNU/Linux 4.4.0-104-generic x86_64)

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/advantage

24 packages can be updated.
0 updates are security updates.


Last login: Tue Jan  9 16:30:24 2018 from 210.75.224.236
[yongxin\@biocloud:~]\$
```

上图中最后一位为命令提示行：中括号内包括`用户名@服务器名：当前所在目录`，`~`代表每个人自己的家目录，如我的是`/mnt/bai/yongxin`

最后一个操作符`\$`即命令提示符，我们可以其后面输入命令来指挥服务器为我们工作啦！！！


### 系统环境介绍

```
pwd # 显示目录位置

mkdir test # 创建名为test的目录

cd test # 进入test目录
```

### 基本文件操作

```
cp /var/www/html/train/usearch/L* ./ # 复制指定目录中以L字母开头的文件至当前目录

ls # 显示当前目录文件

tree # 显示目录结果
```

显示了文件和目录树形结构如下：

```
.
├── L1_1.fq.gz
├── L1_2.fq.gz
├── L1.txt
├── L2_1.fq.gz
├── L2_2.fq.gz
└── L2.txt
```

```
less L1_2.fq.gz # 查看fastq文件右端内容，按q退出预览状态

mv L1_2.fq.gz Library_2.fq.gz

rm L1_1.fq.gz # 删除单个文件

rm L* # 删除所有以L开头的文件

```
### 16S分析流程示例

本流程基于两个文库的fastq数据和相应的实验设计，实现数据分析、Alpha多样性、Beta多样性以及Taxonomy柱状图比较。

```
train_16s.sh

tree
```

分析中所有文件如下：

```
.
├── alpha_boxplot.r
├── alpha_rare_usearch.r
├── beta_cca.r
├── beta_pcoa.r
├── doc
│   ├── design.txt
│   ├── L1.barcode.fa
│   ├── L1.txt
│   ├── L2.barcode.fa
│   └── L2.txt
├── result
│   ├── alpha_rare_groups.pdf
│   ├── alpha_rare_samples.pdf
│   ├── alpha_rare.txt
│   ├── alpha_richness.pdf
│   ├── alpha.txt
│   ├── beta
│   │   ├── bray_curtis.txt
│   │   └── unifrac.txt
│   ├── beta_cca.pdf
│   ├── beta_pcoa_bray_curtis.pdf
│   ├── otus.fa
│   ├── otus.tree
│   ├── otutab_norm.txt
│   ├── otutab_report.txt
│   ├── otutab.txt
│   ├── sintax.txt
│   ├── tax_genus.txt
│   ├── tax_order.txt
│   ├── tax_phylum.txt
│   ├── tax_stack_phylum.pdf
│   └── tax_stack_phylum_sample.png
├── seq
│   ├── L1_1.fq
│   ├── L1_2.fq
│   ├── L2_1.fq
│   ├── L2_2.fq
│   └── rdp_16s_v16.fa
├── tax_phylum.r
└── temp
    ├── filtered.fa
    ├── merge.fq
    ├── stripped.fq
    ├── uniques.fa
    └── zotus.fa
```

主要结果和图表位于`result`目录中，下面我们学习如何从服务器下载文件。


## 上传下载文件Filezilla

我们不仅要学会对Linux发号施令，帮我们工作。更重要的是上传原始测序数据，和下载分析结果。这里我们需要一个软件`FileZilla`。

Filezilla的官方网站：https://filezilla-project.org/

请根据操作系统，下载需要的相应最新版本，并安装好可。

找不到下载链接，可以下载下我预下载的3.30版备用链接。

[Windows 64 bit](http://bailab.genetics.ac.cn/train/FileZilla_3.30.0_win64-setup_bundled.exe)  
[Mac OS X 10.9 or newer](http://bailab.genetics.ac.cn/train/FileZilla_3.30.0_macosx-x86_setup_bundled.dmg)  
[Windows 32 bit](http://bailab.genetics.ac.cn/train/FileZilla_3.30.0_win32-setup_bundled.exe)  

安装成功后运行程序，先点文件菜单下方的`站点管理器`按扭，新建站点，配置服务器信息，连接即可。

(ref:pre-filezilla) 连接到服务器，左侧为本地电脑目录和文件，右侧为服务器目录，鼠标左右键、拖拽操作全支持。

```{r pre-filezilla, fig.cap="(ref:pre-filezilla)"}
knitr::include_graphics("http://bailab.genetics.ac.cn/train/fig/1.03_filezilla.jpg")
```

第一次登陆可能会弹出提醒，勾选“总是信息”，再点确定。即可看到右侧为服务器你的家目录。

在此程序的界面中，不仅可以上传下载文件，还可以对远程和本地文件进行删除、重命名等管理，功能强大。

**下载刚前分析测序数据的结果**

右侧窗口双击进入`test`目录，右键(Mac为双指点击触摸板)选择`result`目录下载即可。注意会动下载至左端当前系统中的目录，想要改变左侧下载位置，请提前选择左侧合适位置。

如果你使用Windows，但无法安装Filezilla，可以选择另一款功能类拟的软件WinSCP 推荐下载[绿色版](http://bailab.genetics.ac.cn/train/WinSCP-5.11.3-Portable.zip)，还可选[安装版WinSCP-5.11.3](http://bailab.genetics.ac.cn/train/WinSCP-5.11.3-Setup.exe)、或访问[主页下载最新版](https://winscp.net/eng/download.php)。


!;
close OUTPUT;



open OUTPUT,">$opts{o}02-usearch.Rmd";
print OUTPUT qq!

# 从原始数据到OTU表 {#usearch}

## 准备测试数据和工作环境 Download test files

此处我准备了两对fastq文件原始数据，和两个txt文件实验设计与每对fastq文件对应。

```
# 1. 准备测试数据和工作环境 Download test files
mkdir test # 创建并进入工作目录
cd test
mkdir -p seq # 原始数据 raw data
mkdir -p doc # 实验设计 design file
mkdir -p temp # 临时文件 temp directory for intermediate files
mkdir -p result # 最终结果 important results
cp /var/www/html/train/usearch/*.gz seq/ # 两个文库的原始测序数据
gunzip seq/* # unzip all file
cp /var/www/html/train/usearch/*.txt doc/ # 实验设计
tree # 查看目录结构和文件 show directories and files in tree format
```

分析起始文件列表如下：

```
── doc
   ├── L1.txt
   └── L2.txt
── result
── seq
   ├── L1_1.fq
   ├── L1_2.fq
   ├── L2_1.fq
   ├── L2_2.fq
   └── rdp_16s_v16.fa
── temp
```

数据分析至少要有两类文件：1类是数据文件，高通量测序数据通常以`.fastq`或`.fq`结尾；另一类是实验设计，通常为文本文件`*.txt`，最重要的是分组信息，对应样品至组，才可以进行组间比较。

Fastq格式解读：FASTQ是基于文本的，保存生物序列（通常是核酸序列）和其测序质量信息的标准格式。其序列以及质量信息都是使用一个ASCII字符标示，最初由Sanger开发，目的是将FASTA序列与质量数据放到一起，目前已经成为高通量测序结果的事实标准。

```
less -S seq/L1_1.fq # 查看fastq文件，示例如下：

\@HISEQ:549:HLYNYBCXY:1:1101:1267:2220 1:N:0:CACTCAAT
TCGTCGCTCGAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGTTGGGCGCTAGGTGTGGGGGACATTCACGTTCTCCGTGCCGTAGCTAACGCATTAAGCGCCCCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAGGAATTGACGGGGACCCGCGCAAGCGGTGGAGCATGTGGTTTAATTCGATGCA
+
DDDDDIHHHIIIIIIIIHIIHIIIIIIIIIIHIHIHIIIIIIIIIIIIIIIIIIIIIIIIIGIIHIHDHHIIHIGHIIIHHHIIIIIII=CHHIEHGCHIIHIIHHIIIH<EHIHGHGHHGHHHIIIIIIIGHIHHIIIIIIII?EHIIIIIGFHIHIC--\@CH<DHCEHEHHEEEHHE?EH/BGH.AGHH?FHHHCEG
\@HISEQ:549:HLYNYBCXY:1:1101:1887:2204 1:N:0:CACTCAAT
TACGAGTATGAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCTACTAGTTGTTGGGTCTTAATTGACTTAGTAACGCAGCTAACGCGTGAAGTAGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAAC
+
DDDD\@H<GHIIIIIIIIIIIIIIIIIIIIHIIIHIIIIIIIIIIGIIIIIIIIFHIHHIIHIIIIIHHIIIFIHIIIIIH?HHHIIIIIIGHIIIIIIHHIDGEHIGIIIHIIIIIIGGGHHHHIGHIIHIIIHHIIHGEGIIIIIIIIIHHIHIIHHCHHHGHCHIIHEHHDD..FHFHHHH.8FHHGEHFHIIIFHH
```

FASTQ文件中每个序列通常有四行：  
序列标识以及相关的描述信息，以‘\@’开头；  
第二行是序列  
第三行以‘+’开头，后面是序列标示符、描述信息，或者什么也不加  
第四行，是质量信息，和第二行的序列相对应，每一个序列都有一个质量评分，根据评分体系的不同，每个字符的含义表示的数字也不相同。  
详见：http://mp.weixin.qq.com/s/1-BjM-hTDlEun2t94On0qg


## 合并双端序列与样品拆分 Merge paired reads and label samples

为什么要双端序列合并，和拆分样品，要从扩增子测序的实验设计説起：

(ref:pip-strcture) 扩增子测序基本结构。Barcode：样品标签，用于混池测序后区分序列来自那个样本；
Primer：在16S/ITS/18S保守区设计的引物，用于扩增rDNA的部分高变区；
Amplicon:扩增的部分 rDNA。

```{r pip-strcture, fig.cap="(ref:pip-strcture)"}
knitr::include_graphics("http://bailab.genetics.ac.cn/train/fig/2.01_amplicon_structure.jpg")
```

(ref:pip-pipe) 测序数据还原为纯净扩增子序列的过程。

```{r pip-pipe, fig.cap="(ref:pip-pipe)"}
knitr::include_graphics("http://bailab.genetics.ac.cn/train/fig/2.02_amplicon_pipe.jpg")
```

### 以L1文库为例实现合并和重命名

```
# 双端序列合并：输入文件，输出文件，重命名序列标签，"\"为了实现代码行，方便阅读
usearch10 -fastq_mergepairs seq/L1_1.fq -reverse seq/L1_2.fq \
	-fastqout temp/L1_merge.fq \
	-relabel L1 # merge, 85.56%

# 质量评估(可选) quality access
usearch10 -fastq_chars temp/L1_merge.fq \
	-log temp/L1_merge.log

# 按barcode标记样品(如果样品已拆分可跳过此步，且mergepairs步加relable sample) label samples

# 提取10 bp的barcode序列
usearch10 -fastx_truncate temp/L1_merge.fq \
	-fastqout temp/L1_barcode.fq \
	-trunclen 10

# 准备barcode样品对应文件 samples barcode in fasta
tail -n+2 doc/L1.txt | cut -f 1-2 | sed 's/^/>/;s/\t/\n/' > doc/L1.barcode.fa 

# 按barcodes与实验设计对应样品重命名序列 label samples
usearch10 -fastx_demux temp/L1_merge.fq -index temp/L1_barcode.fq -barcodes doc/L1.barcode.fa \
	-fastqout temp/L1_demux.fq
```




```
# 2. 合并双端序列与样品拆分 Merge paired reads and label samples
# loop merge and label samples for each library
for i in L1 L2; do
	# Merge paired reads and quality access
	usearch10 -fastq_mergepairs seq/\${i}_1.fq -reverse seq/\${i}_2.fq -fastqout temp/\${i}_merge.fq -relabel \${i} # merge, 85.56%
	usearch10 -fastq_chars temp/\${i}_merge.fq -log temp/\${i}_merge.log # quality access
	# 按barcode标记样品(如果样品已拆分可跳过此步，且mergepairs步加relable \$sample) label samples
	tail -n+2 doc/\${i}.txt | cut -f 1-2 | sed 's/^/>/;s/\\t/\\n/' > doc/\${i}.barcode.fa # 准备barcode样品对应文件 samples barcode in fasta
	usearch10 -fastx_truncate temp/\${i}_merge.fq -trunclen 10 -fastqout temp/\${i}_barcode.fq
	usearch10 -fastx_demux temp/\${i}_merge.fq -index temp/\${i}_barcode.fq -barcodes doc/\${i}.barcode.fa -fastqout temp/\${i}_demux.fq # label samples
done
# 合并多个测序文库 Merge multiple libraries
cat temp/*_demux.fq > temp/merge.fq
```



## 切除引物与质控 Cut primers and quality filter

```
# 3. 切除引物与质控 Cut primers and quality filter
# Cut barcode 10bp + V5 19bp in left and V7 18bp in right
usearch10 -fastx_truncate temp/merge.fq -stripleft 29 -stripright 18 -fastqout temp/stripped.fq
# fastq filter, keep reads error rates less than 1%
usearch10 -fastq_filter temp/stripped.fq -fastq_maxee_rate 0.01 -fastaout temp/filtered.fa
```

## 去冗余与生成OTUs Dereplication and cluster otus

```
# 4. 去冗余与生成OTUs Dereplication and cluster otus
# 去冗余Find unique read sequences and abundances
usearch10 -fastx_uniques temp/filtered.fa -sizeout -relabel U -fastaout temp/uniques.fa
# 预测生物学序列OTU并去除嵌合 Denoise: predict biological sequences and filter chimeras
usearch10 -unoise3 temp/uniques.fa -zotus temp/zotus.fa # 2307 OTUs
sed 's/Zotu/OTU_/g' temp/zotus.fa > result/otus.fa # format OTU prefix
```

## 生成OTU表 Creat OTUs table

```
# 5. 生成OTU表 Creat OTUs table
usearch10 -otutab temp/stripped.fq -otus result/otus.fa -otutabout result/otutab.txt -threads 30 # create OTUs table, 86.2% matched
usearch10 -otutab_stats result/otutab.txt  -output result/otutab_report.txt # Summary OTUs table
usearch10 -otutab_norm result/otutab.txt -sample_size 10000 -output result/otutab_norm.txt # normlize by subsample to 10000
```

我们使用如下命令查看一下OTU表的内容，即获得了每个样品在每种OTUs出现的频次，如果数量少我们是不可手动统计。

```
less -S result/otutab.txt
```

```{r otu-table}
table_design = read.table("result/otutab.txt", sep="\\t", header=T)
knitr::kable(head(table_design), caption="展示OTU表前6行示例。Top 6 lines of OTU table.", booktabs=TRUE)
```


##  Alpha多样性 Alpha diversity

```
# 6. Alpha多样性 Alpha diversity
# Calculate all alpha diversity, details in http://www.drive5.com/usearch/manual/alpha_metrics.html
usearch10 -alpha_div result/otutab_norm.txt -output result/alpha.txt 
# Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch10 -alpha_div_rare result/otutab_norm.txt -output result/alpha_rare.txt  -method without_replacement # 取1%-100%的序列中OTUs数量
```

## Beta多样性 Beta diversity

```
# 7. Beta多样性 Beta diversity
# 基于OTU构建进化树 Make OTU tree
usearch10 -cluster_agg result/otus.fa -treeout result/otus.tree
mkdir result/beta/ # 结果有多个文件，需要目录
# 生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
usearch10 -beta_div result/otutab_norm.txt -tree result/otus.tree -filename_prefix result/beta/
```

## 物种注释 Assign taxonomy

```
# 8. 物种注释 Assign taxonomy
usearch10 -sintax result/otus.fa -db seq/rdp_16s_v16.fa -strand both -tabbedout result/sintax.txt -sintax_cutoff 0.8
# Taxonomy summary reports in phylum, order and genus
usearch10 -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt -rank p -output result/tax_phylum.txt
usearch10 -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt -rank o -output result/tax_order.txt
usearch10 -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt -rank g -output result/tax_genus.txt
# 去除OTUs表头注释符，为R正常识别
sed -i 's/#OTU //' result/otutab*
# Taxonomy中异常字符
sed -i 's/(//g;s/)//g;s/\"//g' result/tax_*.txt
# 合并每个文库实验设计为总表
cat doc/L1.txt <(tail -n+2 doc/L2.txt) > doc/design.txt
```

!;



open OUTPUT,">$opts{o}03-Rplot.Rmd";
print OUTPUT qq!

# 统计绘图

统计绘图，最重要的是在己知实验设计基础上，对各实验组观测数据进行统计与可视化，方便人类观察实验组内、组间相同与差异。

```{r design}
table_design = read.table("doc/design.txt", sep="\\t", header=T)
knitr::kable(table_design, caption="实验设计包括各样品名称、Barcode序列、正反向引物、分组、基本描述等信息，常见的分组还有疾病/健康状态status、批次batch、位置compartment/location/site、时间time/days。Materials information.", booktabs=TRUE)
```

## Alpha多样性指数 Index

```{r alpha}
table_design = read.table("result/alpha.txt", sep="\\t", header=T)
knitr::kable(table_design, caption="各样品常用14种Alpha多样性指数表。Alpha index of each samples.", booktabs=TRUE)
```

以所有样品的Alpha多样性richness(observed OTUs)指数代表物种的多样性，来分析各组间是否存在物种多样性差异。常用的指数还有chao1、dominance和shannon_e。更多指数的详细介绍，请访问 http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html

基于上表中richness列绘制箱线图并统计是否存在差异

```
Rscript alpha_boxplot.r # 基于上表中richness列绘制箱线图并统计是否存在差异
```

```{r alpha_boxplot.r}
# 绘制Alpha稀释曲线

library("agricolae", warn.conflicts = F, quietly = T)
library("ggplot2", warn.conflicts = F, quietly = T)
library("dplyr", warn.conflicts = F, quietly = T)

# 读取实验设计
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

# 读取usearch rarefraction文件
alpha = read.table("result/alpha.txt", header=T, row.names= 1, sep="\t") 

# 提取样品组信息
sampFile = as.data.frame(design\$genotype,row.names = row.names(design))
colnames(sampFile)[1] = "group"

## richness index
# add design to alpha
index = cbind(alpha[rownames(design),]\$richness, sampFile) 
colnames(index) = c("richness","group") # add richness colname is value
# 统计各组间差异
model = aov(richness ~ group, data=index)
out <- LSD.test(model,"group", p.adj="none")
stat = out\$groups
# 分组结果添入Index
index\$stat=stat[as.character(index\$group),]\$groups
# 设置分组位置为各组y最大值+高的3%
max=max(index\$richness)
min=min(index\$richness)
x = index[,c("group","richness")]
y = x%>% 
  group_by(group)%>%
  summarise(Max=max(richness))
y=as.data.frame(y)
rownames(y)=y\$group
index\$y=y[as.character(index\$group),]\$Max + (max-min)*0.03

p = ggplot(index, aes(x=group, y=richness, color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y="richness index") + theme_classic() +
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
p
```

 箱线图展示各样品及组的Alpha多样性分布，方法采用 Richness (Observed OTUs) index，只包括样品OTU种类信息。图中KO(knock out)代表基因敲除突变体，OE(over expression)代表基因过表达株系，WT(wild type)代表野生型。各组间采用R语言agricolae包的LSD.test函数统计，两组上方无相同字母代表组间存在显著差异(Pvalue < 0.05)。

## Alpha稀释曲线 Rarefraction curve

```{r Rarefraction}
table_design = read.table("result/alpha_rare.txt", sep="\\t", header=T)
knitr::kable(head(table_design), caption="OTU表按百分比重抽样后统计丰富度表前6行。Top 6 lines of richness index in rarefraction curve.", booktabs=TRUE)
```

运行绘图的脚本：

```
Rscript alpha_rare_usearch.r 
```

```{r alpha_rare_usearch.r }
# 绘制Alpha稀释曲线

library("reshape2", warn.conflicts = F, quietly = T)
library("ggplot2", warn.conflicts = F, quietly = T)

# 读取实验设计
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

# 读取usearch rarefraction文件
rare = read.table("result/alpha_rare.txt", header=T, row.names= 1, sep="\t") 
# 提取样品组信息
sampFile = as.data.frame(design\$genotype,row.names = row.names(design))
colnames(sampFile)[1] = "group"

# 直接展示样品
rare\$x = rownames(rare) # 添加x轴列
rare_melt = melt(rare, id.vars=c("x")) # 转换为长表格
rare_melt\$x = factor(rare_melt\$x, levels=1:100) # 设置x轴顺序

rare_melt3 = merge(sampFile,rare_melt, by.x="row.names", by.y="variable")
rare_melt3\$variable=rare_melt3\$Row.names

# 按样品分组，按组上色
p = ggplot(rare_melt3, aes(x = x, y = value, group = variable, color = group )) + 
  geom_line()+xlab("Rarefraction Percentage")+ylab("Richness (Observed OTUs)")+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+ theme_classic()
p
ggsave(paste("result/alpha_rare_samples.pdf", sep=""), p, width = 8, height = 5)
ggsave(paste("result/alpha_rare_samples.png", sep=""), p, width = 8, height = 5)


# 求各组均值
# 读取usearch rarefraction文件，上面己经修改，必须重新读入
rare = read.table("result/alpha_rare.txt", header=T, row.names= 1, sep="\t") 
# 文件的行名为纯数字，转置将整个矩阵num变为字符chr，需要提前修改行名；是因为上方修改了格式，重读解决
# rownames(rare)=c(paste("a",1:100,sep=""))
# 转置rare表格与实验设计合并，并去除第一列样品名
mat_t = merge(sampFile, t(rare), by="row.names")[,-1]
# 按第一列合并求均值
mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
# 修正行名
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno

rare=as.data.frame(round(mat_mean_final))
rare\$x = rownames(rare)
rare_melt = melt(rare, id.vars=c("x"))
#rare_melt\$x = factor(rare_melt\$x, levels=1:100)

# 求各组标准误
# 转置rare表格与实验设计合并，并去除第一列样品名
se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se) # se 为什么全是NA
mat_se_final = do.call(rbind, mat_se)[-1,]
colnames(mat_se_final) = geno

rare_se=as.data.frame(round(mat_se_final))
rare_se\$x = rownames(rare_se)
rare_se_melt = melt(rare_se, id.vars=c("x"))

# 添加标准误到均值中se列
rare_melt\$se=rare_se_melt\$value
# 去除之前添加的x轴字符，转换为纯数字，第一步样品绘图修改了原始数据，第二步要重新读入
# rare_melt\$x1=gsub("a","",rare_melt ,perl=TRUE)
# 添加levels顺序，否则
rare_melt\$x = factor(rare_melt\$x, levels=c(1:100))

p = ggplot(rare_melt, aes(x = x, y = value, group = variable, color = variable )) + 
  geom_line()+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
  xlab("Percentage")+ylab("Richness (Observed OTUs)")+theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=1, hjust=1))+
  scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10) 
p
ggsave(paste("result/alpha_rare_groups.pdf", sep=""), p, width = 8, height = 5)
ggsave(paste("result/alpha_rare_groups.png", sep=""), p, width = 8, height = 5)
```

ggplot2可视化样品(按组着色)和组(组内样品求均值)稀释曲线图；Rarefraction curve of richness in samples and groups.

## Beta主坐标轴分析 PCoA

```
Rscript beta_pcoa.r
```

```{r beta_pcoa.r}
# 绘制主坐标轴分析PCoA

library("vegan", warn.conflicts = F, quietly = T)
library("ggplot2", warn.conflicts = F, quietly = T)

# 读取实验设计
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

#  PCoA bray_curtis
bray_curtis = read.table("result/beta/bray_curtis.txt", sep="\t", header=T,  row.names= 1)

# subset matrix and design
# idx = rownames(sub_design) %in% colnames(bray_curtis) 
# sub_design = sub_design[idx,]
# bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design)] # subset and reorder distance matrix

# cmdscale {stats}, Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa\$points) # get coordinate string, format to dataframme
eig = pcoa\$eig
points = cbind(points, design[rownames(points),]\$genotype)
colnames(points) = c("x", "y", "z","group") 

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group)) + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="bray_curtis PCoA")  + stat_ellipse(level=0.68) + theme_classic()
p
ggsave("result/beta_pcoa_bray_curtis.pdf", p, width = 8, height = 5)
ggsave("result/beta_pcoa_bray_curtis.png", p, width = 8, height = 5)

```

## 限制性主坐标轴分析 Constrained PCoA

```
Rscript beta_cca.r
```

```{r beta_cca.r}
library("vegan", warn.conflicts = F, quietly = T)
library("ggplot2", warn.conflicts = F, quietly = T)

# 读取实验设计
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 
design\$group=design\$genotype
#  PCoA bray_curtis
sub_otu_table = read.table("result/otutab_norm.txt", sep="\t", header=T,  row.names= 1)

variability_table = function(cca){
  chi = c(cca\$tot.chi, cca\$CCA\$tot.chi, cca\$CA\$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}


# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(sub_otu_table) ~ group, data=design, add=F, sqrt.dist=T, distance="bray") 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen, permutations = 10000, parallel = 9)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen\$CCA\$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen\$CCA\$wa[, 1:2]
points = as.data.frame(points)
colnames(points) = c("x", "y")
points = cbind(points, design[match(rownames(points), rownames(design)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group, shape=group))+
  geom_point() +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + 
  theme_classic()+ stat_ellipse(level=0.68)
p
```

## 物种丰度柱状图 Taxonomy barplot

```
Rscript tax_phylum.r
```

```{r tax_phylum.r}
library("reshape2", warn.conflicts = F, quietly = T)
library("ggplot2", warn.conflicts = F, quietly = T)

# 读取实验设计
design = read.table("doc/design.txt", header=T, row.names= 1, sep="\t") 

#  PCoA bray_curtis
sample_tax = read.table("result/tax_phylum.txt", sep="\t", header=T,  row.names= 1)
sample_tax=sample_tax[,rownames(design)]

# Stackplot for each samples
#sample_tax = read.delim("sum_taxa/otu_table_tax_L2.txt", row.names= 1,  header=T, sep="\t")
mean_sort = sample_tax[(order(-rowSums(sample_tax))), ] # decrease sort
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")
if (TRUE){
  tax_name = gsub("[\\\\w;_]+__","",rownames(mean_sort),perl=TRUE)
  j=1
  for (i in 1:length(tax_name)){
    if (tax_name[i]==""){
      tax_name[i]=paste("Noname",j,sep='')
      j=j+1
    }
  }	
  rownames(mean_sort) = tax_name # rowname unallowed same name
}
merge_tax=mean_sort
mean_sort\$phylum = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("phylum")))
# data_all\$phylum  = factor(data_all\$phylum, levels=rownames(mean_sort))   # set taxonomy order
data_all = merge(data_all, design[c("genotype")], by.x="variable", by.y = "row.names")

p = ggplot(data_all, aes(x=variable, y = value, fill = phylum )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ genotype, scales = "free_x", switch = "x") +  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p

# 按组均值绘制柱状图

sampFile = as.data.frame(design\$genotype,row.names = row.names(design))
colnames(sampFile)[1] = "group"
mat_t = t(merge_tax)

mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean\$group
colnames(mat_mean_final) = geno

mean_sort=as.data.frame(mat_mean_final)
mean_sort\$phylum = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("phylum")))
# data_all\$phylum  = factor(data_all\$phylum, levels=rownames(mean_sort))   # set taxonomy order
# data_all = merge(data_all, design[c("genotype")], by.x="variable", by.y = "row.names")

p = ggplot(data_all, aes(x=variable, y = value, fill = phylum )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()
p
```

!;

close OUTPUT;

#- 课题目的Aim
#
#$list{aim}
#
#- 材料描述Materials
#
#样品材料描述信息详见 Table \\\@ref(tab:material) 所示，[完整实验设计表material.txt](doc/material.txt)。
#
#```{r material}
#table_design = read.table("doc/material.txt", sep="\\t", header=T)
#knitr::kable(table_design, caption="各材料/样品组的原始名、新名称、基本描述等信息。Materials information.", booktabs=TRUE)
#```
#
#- 样品列表Mapping file
#
#$list{design}样品信息如 Table \\\@ref(tab:design) 所示，[完整实验设计表design.xls](doc/design.xls)。
#
#```{r design}
#table_design = read.table("doc/design.txt", sep="\\t", header=T)
#table_design=head(table_design,n=10)
#knitr::kable(table_design, caption="样品信息前10行。Top 10 lines of sample information.", booktabs=TRUE)
#```


open OUTPUT,">$opts{o}98-scheme.Rmd";
print OUTPUT qq!
# 方案Scheme {#project_scheme}

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



#if (-e $opts{l}) {
#
### 读取文库列表文件
#open DATABASE,"<$opts{l}";
#<DATABASE>;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @library,$tmp[0] if defined($tmp[0]); # 判断是否非空的行
#}
#close DATABASE;
#
#open OUTPUT,">$opts{o}04-a-sequenceQuality.Rmd";
#print OUTPUT qq!
## 质量控制Quality control {#sequencing_quality_summary}
#
### 质量评估Quality assessment {#sub-sequence-qc}
#
#说明：16S扩增子测序数据主要来自HiSeq2500产出的双端各250 bp (PE250)数据，优势是读长长且价格适中(性价比高)。HiSeqX PE150和MiSeq PE300也比较常见，但PE150过短分辨率低，而PE300价格高且末端序列质量过低。此外454在之前研究较多且设备已经停产，PacBio读长长可直接测序16S全长1.5kb代表未来的趋势。测序公司通常会返回raw data和clean data两种数据，raw data为测序获得的原始数据，而clean data则为去除含有接头序列及测序不确定N比例较高的结果，通常直接采用clean data进行质量评估及后续分析。数据质量评估结果中测序reads碱基质量分布图，常用于判断建库或测序有无异常。序列重复情况分布，判断原始序列的DNA质量、重复序列比例及PCR扩增重复情况，如重复序列较高可能某些菌高丰度或PCR扩增导致，对低丰度菌的结果影响较大。
#
#!;
#close OUTPUT;
#
#foreach $library (@library) {
#open OUTPUT,">>$opts{o}04-a-sequenceQuality.Rmd";
##`cp -r clean_data/${library}_*.html result/`;
##`cp -r clean_data/${library}_1_fastqc/ result/`;
##`cp -r clean_data/${library}_2_fastqc/ result/`;
##(ref:quality-fastqc-${library}) 测序Reads质量评估文库${library}。Clean reads左端(A)和右端(B)数据质量评估；clean reads左端(C)和右端(D)序列重复情况分布。Quality control of clean reads [HTML report of library ${library}_1](result/${library}_1_fastqc.html)  [HTML report of library ${library}_2](result/${library}_2_fastqc.html)
##```{r quality-fastqc-${library}, fig.cap="(ref:quality-fastqc-${library})", out.width="49%"}
##figs_1 = paste0("result/${library}_", c("1_fastqc/Images/per_base_quality", "2_fastqc/Images/per_base_quality", "1_fastqc/Images/duplication_levels", "2_fastqc/Images/duplication_levels"),".png")
##knitr::include_graphics(figs_1)
##```
#print OUTPUT qq!
#**文库${library}质量评估Quality assessment**
#
#文库${library}测序数据质量评估报告。Quality assessment clean reads of library ${library}。
#
#- [Illumina测序左端数据质量评估报告，HTML report of library ${library}_1](clean_data/${library}_1_fastqc.html)  
#- [Illumina测序右端数据质量评估报告，HTML report of library ${library}_2](clean_data/${library}_2_fastqc.html)
#
#(ref:quality-split-${library}) 文库${library}各样品按barcode拆分获得的高质量测序数据，按实验设计组着色。Distribution of sequenced reads of samples in library ${library}. Samples were colored by group information. 1 Million = 10^6^. [PDF](result/stat_lib_split_${library}.pdf)
#
#```{r quality-split-${library}, fig.cap="(ref:quality-split-${library})"}
#knitr::include_graphics("result/stat_lib_split_${library}.png")
#```
#
#!;
#close OUTPUT;
#}
#
#open OUTPUT,">>$opts{o}04-a-sequenceQuality.Rmd";
#`ln result/qc.sum result/qc.xls -f`;
#`echo -e "Library\tCount\tLength">temp/temp`;
#`cat temp/temp result/length.txt > result/length.xls`;
#`rm temp/temp`;
#print OUTPUT qq!
### 质控流程Quality control flow {#sub-sequence-summary}
#
#(ref:quality-sum) 测序文库数据量和长度分布。(A) 柱状图展示各文库数据标准化筛选各步骤有效数据分布。主要包括数据低质量及污染序列过滤、双端合并、筛选扩增子并统一序列方向、按barcode拆分样品、去除5’引物序列、去除3’引物序列为下一步分析的高质量样本序列；(B) 折线图展示各测序文库中序列的长度分析。Data size and length distribution of sequencing libraries. (A) Bar plot showing reads count of each library in read filter process. (B) Line plot showing reads length distribution of each library. [Sum PDF](result/stat_lib_qc_sum.pdf) [Sum XLS](result/qc.xls) [Length PDF](result/stat_lib_length.pdf) [Length XLS](result/length.xls)
#
#```{r quality-sum, fig.cap="(ref:quality-sum)"}
#knitr::include_graphics(c("result/stat_lib_qc_sum.png","result/stat_lib_length.png"))
#```
#
#!;
#close OUTPUT;
#
#}
#
#open OUTPUT,">$opts{o}04-b-tree.Rmd";
#print OUTPUT qq!
## OTU概述Summary {#result-tree}
#
### OTU概述Summary  {#sub-result-table}
#
#```{r otutable}
#tab = read.table("result_k1-c/otu_table.txt", sep="\\t", header=T)
#tab=head(tab,n=10)
#knitr::kable(tab, caption="OTU表前10行。Top 10 lines of OTU table.", booktabs=TRUE)
#```
#
#原始数据经过去冗余、只保留出现8次以上的序列、使用unoise3进行去噪(非97%聚类，相当于100%聚类)、自身比对去嵌合体、与RDP数据比对去嵌合体、与GreenGene 13.5(可选RDP11.5/SILVA128)比对去除非细菌序列。原始含分类学注释OTU表下载[BIOM](result/otu_table_tax.biom) [TXT](result/otu_table_tax.txt)，[代表性序列Fasta](result/rep_seqs.fa)。基本统计信息如下：
#
#```{r otutable-sum}
#tab = read.table("result/otu_table.sum", sep=":", header=F)
#colnames(tab)=c("Types","Values")
#tab=head(tab,n=10)
#knitr::kable(tab, caption="OTU表基本统计信息。Summary of OTU table.", booktabs=TRUE)
#```
#
#原始OTU表去除Cyanobacteria和Chloroflexi菌门细菌(可能为叶绿体/线粒体来源，ITS去除宿主)后，只保留在至少在一个样品中存在丰度大于$opts{a}的OTU，接下来的数据进行差异统计分析。Raw OTU remove Cyanobacteria and  abundance lower than $opts{a}. [BIOM](result_k1-c/otu_table_tax.biom) [TXT](result_k1-c/otu_table_tax.txt) [Normalized OTU table in percentage](result_k1-c/otu_table_norm.txt)，[代表性序列Fasta](result_k1-c/rep_seqs.fa)。Basic statistics information：
#
#```{r otutable-filter-sum}
#tab = read.table("result_k1-c/otu_table.sum", sep=":", header=F)
#colnames(tab)=c("Types","Values")
#tab=head(tab,n=10)
#knitr::kable(tab, caption="OTU表基本统计信息。Summary of OTU table.", booktabs=TRUE)
#```
#
### 物种分类树Taxonomy {#sub-result-graphlan}
#
#(ref:tree-graphlan) 高丰度OTU(>0.5%)物种注释分类树，按分类学目(A. order)、科(B. family)、属(C. genus)水平进行文本标签注释，结果可以看到本实验中鉴定的细菌OTU主要分布于不同分类级别的哪些目、科、属。Taxonomy tress show high abundance OTU (>0.5%), and their order, family and genus annotaion of taxonomy. [order PDF](result/tax_order.pdf)  [family PDF](result/tax_family.pdf)  [genus PDF](result/tax_genus.pdf)
#
#```{r tree-graphlan, fig.cap="(ref:tree-graphlan)"}
#if ($opts{S}) {figs_2 = paste0("result/tax_", c("order"),".png")
#}else{figs_2 = paste0("result/tax_", c("order", "family", "genus"),".png")}
#knitr::include_graphics(figs_2)
#```
#
### 进化树Phylogenetic {#sub-result-ggtree}
#
#(ref:tree-ggtree) 高丰度OTU系统发生树分析(>0.5%)，按分类学门(A. phylum)、纲(B. class)、目(C. order)水平进行着色(不同分类图请点击段末链接)，结果可以看到本实验中鉴定的细菌OTU主要分布于那些分类级别，同时表明细菌的物种注释信息与16S的序列发生树的进化关系高度一致(图中显示的0组为存为分类不明确的组，请使用时用AI删除相应图注)。Phylogenetic tress show high abundance OTU (>0.5%), and their phylum, class and order annotaion of taxonomy. [phylum PDF](result/ggtree_phylum.pdf); [class PDF](result/ggtree_class.pdf); [order PDF](result/ggtree_order.pdf).
#
#```{r tree-ggtree, fig.cap="(ref:tree-ggtree)", out.width="99%"}
#if ($opts{S}) {figs_1 = paste0("result/ggtree_", c("phylum"),".png")
#}else{figs_1 = paste0("result/ggtree_", c("phylum", "class", "order"),".png")}
#knitr::include_graphics(figs_1)
#```
#!;
#close OUTPUT;
#
#
#open OUTPUT,">$opts{o}04-d-diversity.Rmd";
#`ln result/alpha.txt result/alpha.xls -f`;
#`ln result_k1-c/beta.txt result_k1-c/beta.xls -f`;
#print OUTPUT qq!
## 多样性Diversity {#result-diversity}
#
### 相关分析Correlation {#result-diversity-cor}
#
#(ref:group-cor) 基于各样品组OTU均值计算皮尔森相关系数。Pearson correlation of all groups (mean). [PDF](result_k1-c/heat_cor_groups.pdf)
#
#```{r group-cor, fig.cap="(ref:group-cor)", out.width="99%"}
#knitr::include_graphics("result_k1-c/heat_cor_groups.png")
#```
#
#(ref:sample-cor) 基于各样品组OTU相对丰度计算皮尔森相关系数。Pearson correlation of all samples. [PDF](result_k1-c/heat_cor_samples.pdf)
#
#```{r sample-cor, fig.cap="(ref:sample-cor)", out.width="99%"}
#knitr::include_graphics("result_k1-c/heat_cor_samples.png")
#```
#
### Alpha多样性 {#result-diversity-alpha}
#
#主要展示各样品、组的物种丰富度(richness)、均匀度(evenness)的分布
#
#### Alpha diversity value 多样性数值 {#result-diversity-alpha-value}
#
#各样品常用四种Alpha多样性计算方法结果见如 Table \\\@ref(tab:alpha) 所示。[完整表格下载TXT](result/alpha.xls)
#
#```{r alpha}
#table_alpha = read.table("result/alpha.xls", sep="\\t", header=T, row.names = 1)
#table_alpha=head(table_alpha,n=10)
#knitr::kable(table_alpha, caption="样品四种Alpha多样性结果前10行", booktabs=TRUE)
#```
#
#说明：常见的Alpha多样性相关介绍和背景知识，推荐阅读[《箱线图：Alpha多样性，老板再也不操心我的文献阅读》](https://mp.weixin.qq.com/s/CkHVLzDVosKzoxFEIUpNWw)一文，更多相关知识可阅读[Alpha diversity measures](http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html)。  
#
#
#### Alpha多样性分布图 {#result-diversity-alpha-figure}
#
#[QIIME绘制稀释曲线 Rarefraction plots](result/rarefaction_plots.html)
#
#(ref:div-alpha) 箱线图展示各样品及组的微生物组Alpha多样性，方法采用(A) Shannon index，包括样品的可操作分类单元(operational taxonomic unit, OTU)种类(richness)及丰度(evenness)信息；(B) Observed OTUs index，只包括样品OTU种类信息。(C) Chao1 index,基于样品测序中单拷贝OTU(饱合情况)估算样品物种种类的方法; (D) PD whole tree index, 多样性评估时考虑OTU间的进化关系，通常进化关系相近的物种可能存在丰度更相关。图中KO(knock out)代表基因敲除突变体，OE(overexpression)代表基因过表达株系，WT(wild-type)代表野生型。各组间分组采用R语言agricolae包的LSD.test函数统计，两组上方无相同字母代表组间存在显著差异(pvalue < 0.05)。附文本有t-test方法统计各组间是否存在显著差异的p-value水平。
#[Shannon TXT](result/alpha_shannon_stats.txt)  [observed_otus TXT](result/alpha_observed_otus_stats.txt)  [chao1 TXT](result/alpha_chao1_stats.txt)  [PD_whole_tree TXT](result/alpha_PD_whole_tree_stats.txt)
#Within sample diversity (α-diversity) measurements among each genotype. (A) Shannon index, estimated species richness and evenness; (B) Observed OTUs index, only calculate species richness; (C) Chao1 index, calculate richness based on observed, singletons and doubletons; (D) PD whole tree index, diversity considered the evolution distance as weighted. These results indicate genotype not significantly change microbial diversity. The horizontal bars within boxes represent median. The tops and bottoms of boxes represent 75th and 25th quartiles, respectively. The upper and lower whiskers extend 1.5× the interquartile range from the upper edge and lower edge of the box, respectively. All outliers are plotted as individual points (Edwards et al., 2015).
# [Shannon PDF](result/alpha_shannon.pdf)  [observed_otus PDF](result/alpha_observed_otus.pdf)  [chao1 PDF](result/alpha_chao1.pdf)  [PD_whole_tree PDF](result/alpha_PD_whole_tree.pdf) 
#
#```{r div-alpha, fig.cap="(ref:div-alpha)", out.width="99%"}
##figs_2 = paste0("result/alpha_", c("shannon", "observed_otus", "chao1", "PD_whole_tree"),".png")
#if ($opts{S}) {figs_2 = paste0("result/alpha_", c("shannon", "observed_otus"),".png")
#}else{figs_2 = paste0("result/alpha_", c("shannon", "observed_otus", "chao1", "PD_whole_tree"),".png")}
#knitr::include_graphics(figs_2)
#```
#
#
### PCoA样品间差异Beta diversity
#
#(ref:div-beta) 主坐标轴分析(PCoA)展示第1/2坐标轴下各样品间微生物组差异(dissimilarity)，距离计算方法采用(A) bray curtis; (B) unweighted unifrac; (C) weighted unifrac。[采用Adonis统计各样品组间的显著性差异P值](result_k1-c/beta.xls)。
#Principal coordinate analysis (PCoA) using the (A) bray curtis metric, (B) unweighted unifrac metric and (C) weighted unifrac metric shows dissimilarity of microbial communities. [bray_curtis PDF](result_k1-c/beta_pcoa_bray_curtis.pdf)  [unweighted_unifrac PDF](result_k1-c/beta_pcoa_unweighted_unifrac.pdf)  [weighted_unifrac PDF](result_k1-c/beta_pcoa_weighted_unifrac.pdf)  [bray_curtis with sampleID lables PDF](result_k1-c/beta_pcoa_bray_curtis_lab.pdf)
#
#```{r div-beta, fig.cap="(ref:div-beta)", out.width="99%"}
##figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis", "unweighted_unifrac", "weighted_unifrac"),".png")
#if ($opts{S}) {figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis"),".png")
#}else{figs_2 = paste0("result_k1-c/beta_pcoa_", c("bray_curtis", "unweighted_unifrac", "weighted_unifrac"),".png")}
#knitr::include_graphics(figs_2)
#```
#!;
#
## 检测限制条件主成分分析是否有结果图，没有则不输出
#$file = "result_k1-c/CPCoA_$opts{g}.png";
#if (-e $file) {
#print OUTPUT qq!
#
### 限制性主坐标轴分析(全部OTU)CPCoA (Total OTU) 
#
#(ref:div-CPCoA-a) 以基因型为条件分析其贡献率和样品组间差异。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
#Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result/CPCoA_$opts{g}.pdf)  [PDF labels](result/CPCoA_$opts{g}_lab.pdf)  
#
#```{r div-CPCoA-a, fig.cap="(ref:div-CPCoA-a)", out.width="99%"}
#knitr::include_graphics("result/CPCoA_$opts{g}.png")
#```
#	
### 限制性主坐标轴分析(高丰度OTU)CPCoA (High abundance OTU) 
#
#(ref:div-CPCoA) 以基因型为条件分析其贡献率和样品组间差异(筛选至少在一个样品中OTU丰度 > $opts{a} )。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
#Constrained principal coordinate analysis on bacterial microbiota (Only OTU abundance more than $opts{a} in one sample were kept). Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result_k1-c/CPCoA_$opts{g}.pdf)  [PDF labels](result_k1-c/CPCoA_$opts{g}_lab.pdf)  
#
#
#```{r div-CPCoA, fig.cap="(ref:div-CPCoA)", out.width="99%"}
#knitr::include_graphics("result_k1-c/CPCoA_$opts{g}.png")
#```
#
#!;
#
#
## 检测第二分组是否有多个条件，有则逐个条件输出
#$opts{D}=~s/"//g;
#@g2_list=split(/,/,$opts{D});
#if (@g2_list>1) {
#foreach $i(@g2_list) {
#
#$file = "result_k1-c/CPCoA_$opts{g}_${i}.png";
#if (-e $file) {
#print OUTPUT qq!
#### 限制性主坐标轴分析${i}组
#
#(ref:div-CPCoA-${i}) 在${i}组内，以基因型为条件分析其贡献率和样品组间差异。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
#Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result_k1-c/CPCoA_$opts{g}_${i}.pdf)  
#
#```{r div-CPCoA-${i}, fig.cap="(ref:div-CPCoA-${i})", out.width="99%"}
#knitr::include_graphics("result_k1-c/CPCoA_$opts{g}_${i}.png")
#```
#!;
#}
#}
#}
## 检测第三分组是否有多个条件，有则逐个条件输出
#$opts{F}=~s/"//g;
#@g3_list=split(/,/,$opts{F});
#if (@g3_list>1) {
#foreach $i(@g3_list) {
#
#$file = "result_k1-c/CPCoA_$opts{g}_${i}.png";
#if (-e $file) {
#print OUTPUT qq!
#### 限制性主坐标轴分析${i}组
#
#(ref:div-CPCoA-${i}) 在${i}组内，以基因型为条件分析其贡献率和样品组间差异。vriance代表当前基因型条件下各样品间差异所占的比重或贡献率，P值示基因型各组间是否存在显著差异，各样品间距离计算方法为Bray-Curtis distances。
#Constrained principal coordinate analysis on bacterial microbiota. Variation between samples in Bray-Curtis distances constrained by genotype. (Bulgarelli et al., 2015).[PDF](result_k1-c/CPCoA_$opts{g}_${i}.pdf)  
#
#```{r div-CPCoA-${i}, fig.cap="(ref:div-CPCoA-${i})", out.width="99%"}
#knitr::include_graphics("result_k1-c/CPCoA_$opts{g}_${i}.png")
#```
#!;
#}
#}
#}
#}
#close OUTPUT;
#
#
#
### 样品组比较
#open DATABASE,"<$opts{c}";
#while (<DATABASE>) {
#	chomp;
#	push @group,$_;
#}
#close DATABASE;
#@tax_en=qw#phylum class order family genus#;
#@tax_cn=qw#门 纲 目 科 属#;
#
#open OUTPUT,">$opts{o}04-e-taxonomy.Rmd";
#print OUTPUT qq!
#
## 高分类级别差异Different Taxonomy {#result-taxonomy}
#
### 差异分类学单元数量Different taxonomy summary 
#
#样品组在不同分类级别上显著差异分类单元数量。Different taxonomy unit under each level (Pvalue < 0.05, FDR < 0.05)如 Table \\\@ref(tab:taxonomy-sum) 所示。[TXT](result_k1-c/tax_sum.txt)
#
#```{r taxonomy-sum}
#table_taxonomy = read.table("result_k1-c/tax_sum.txt", sep="\t", header=T)
#knitr::kable(table_taxonomy, caption="样品组显著差异taxonomy数量", booktabs=TRUE)
#```
#
#!;
#
#foreach $i (0..4) {
#
## 开启精简模式，只输入门、目水平差异，纲、科和属跳过
#if ($opts{S} eq "TRUE") {
#	next if $i==1;
#	next if $i==2;
#	next if $i==4;
#}
#
#print OUTPUT qq!
### $tax_cn[$i]水平组间差异$tax_en[$i] {#result-taxonomy-$tax_en[$i]}
#
#(ref:taxonomy-$tax_en[$i]) 柱状图展示各样品组微生物组分类学$tax_cn[$i]水平相对丰度。(A) 堆叠柱状图展示各样品相对丰度。(B) 堆叠柱状图展示各组平均相对丰度，X轴为各样品组，Y轴为各$tax_cn[$i]类相对百分比，只列出了丰度大于0.1%的$tax_cn[$i]，其它所有$tax_cn[$i]归入Low Abundance类。(C) 条形图展示最高丰度的五大菌$tax_cn[$i]平均丰度及标准误，我们可以观察各样品组$tax_cn[$i]水平上相关丰度的差异及组内生物学重复间的波动范围。
#Bar plot showing $tax_en[$i] abundances in each genotype. (A) Stack plot showing high abundance (>0.1%) $tax_en[$i] in each sample;  (B) Stack plot showing high abundance (>0.1%) $tax_en[$i] in each group; (C) Bar plot showing top 5 $tax_en[$i] abundance and error bar in each genotype.  [stack sample PDF](result_k1-c/tax_stack_$tax_en[$i]_sample.pdf)  [stack group PDF](result_k1-c/tax_stack_$tax_en[$i]_top9.pdf)  [bar PDF](result_k1-c/tax_bar_$tax_en[$i]_top5.pdf) [raw Data](result_k1-c/database_$tax_en[$i].txt)
#
#```{r taxonomy-$tax_en[$i], fig.cap="(ref:taxonomy-$tax_en[$i])", out.width="99%"}
##figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_top9", "bar_$tax_en[$i]_top5"),".png")
#if ($opts{S}) {figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_sample","stack_$tax_en[$i]_top9"),".png")
#}else{figs_2 = paste0("result_k1-c/tax_", c("stack_$tax_en[$i]_sample","stack_$tax_en[$i]_top9", "bar_$tax_en[$i]_top5"),".png")}
#knitr::include_graphics(figs_2)
#```
#!;
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#$file = "result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.pdf";
##print -e $file,"\n"; # 不存在为uninitialized，存在返回1
#
#
#if (-e $file) {
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.txt > result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls`;
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.txt > result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls`;
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1]
#
#$tmp[0]与$tmp[1]相比显著差异的分类单元信息如 Table \\\@ref(tab:taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]) 所示。[Enriched TXT](result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls)
#
#```{r taxonomy-$tmp[0]vs$tmp[1]-$tax_en[$i]}
#table_taxonomy_e = read.table("result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
#table_taxonomy_d = read.table("result_k1-c/$tax_en[$i]_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
#table_taxonomy_merge = rbind(table_taxonomy_e,table_taxonomy_d)
##table_taxonomy = table_taxonomy_merge[,1:5]
##table_taxonomy = table_taxonomy_merge
#table_taxonomy = head(table_taxonomy_merge,n=10)
#knitr::kable(table_taxonomy, caption="样品组$tmp[0] vs $tmp[1]显著差异$tax_cn[$i]前10行；Significantlly different $tax_en[$i].", booktabs=TRUE)
#```
#
#(ref:taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1]) 热图展示$tmp[0]vs$tmp[1]在$tax_cn[$i]水平差异分类单元。Heatmap show differentially abundance $tax_en[$i].[PDF](result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.pdf)
#
#```{r taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1], fig.cap="(ref:taxonomy-$tax_en[$i]-$tmp[0]vs$tmp[1])", out.width="99%"}
#knitr::include_graphics("result_k1-c/heat_$tax_en[$i]_$tmp[0]vs$tmp[1]_sig.png")
#```
#
#!;
#}else{
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1]
#
#无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.
#
#!;
#}
#}
#
#}
#
#
#
## phylum + proteobacteria class
#print OUTPUT qq!
### 门及变形菌纲差异 Phylum and class of Proteobacteria  {#result-taxonomy-phylumpro}
#
#(ref:taxonomy-phylumpro) 柱状图展示各样品组微生物组分类学门及变形菌纲水平相对丰度。(A) 堆叠柱状图展示各组平均相对丰度，X轴为各样品组，Y轴为各门及变形菌纲类相对百分比，只列出了丰度大于0.1%的门及变形菌纲，其它所有门及变形菌纲归入Low Abundance类。(B) 叠柱状图展示各样品相对丰度，我们可以观察各组内样品间的波动情况，又可以看到组间差异。
#Bar plot showing phylumpro abundances in each genotype. (A) Stack plot showing high abundance (>0.1%) phylumpro; (B) Stack plot showing high abundance (>0.1%) phylumpro in each sample. [stack PDF](result_k1-c/tax_stack_phylumpro_abc.pdf)  [stack sample PDF](result_k1-c/tax_stack_phylumpro_sample.pdf)  [raw Data](result_k1-c/database_phylumpro.txt)
#
#```{r taxonomy-phylumpro, fig.cap="(ref:taxonomy-phylumpro)", out.width="99%"}
##figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc", "stack_phylumpro_sample"),".png")
#if ($opts{S}) {figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc"),".png")
#}else{figs_2 = paste0("result_k1-c/tax_", c("stack_phylumpro_abc", "stack_phylumpro_sample"),".png")}
#knitr::include_graphics(figs_2)
#```
#!;
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#$file = "result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.pdf";
##print -e $file,"\n"; # 不存在为uninitialized，存在返回1
#
#
#if (-e $file) {
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.txt > result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls`;
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.txt > result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls`;
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1]
#
#$tmp[0]与$tmp[1]相比显著差异的分类单元信息如 Table \\\@ref(tab:taxonomy-$tmp[0]vs$tmp[1]-phylumpro) 所示。[Enriched TXT](result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls)
#
#```{r taxonomy-$tmp[0]vs$tmp[1]-phylumpro}
#table_taxonomy_e = read.table("result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
#table_taxonomy_d = read.table("result_k1-c/phylumpro_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
#table_taxonomy_merge = rbind(table_taxonomy_e,table_taxonomy_d)
#table_taxonomy = head(table_taxonomy_merge,n=10)
#knitr::kable(table_taxonomy, caption="样品组$tmp[0] vs $tmp[1]显著差异门及变形菌纲前10行；Significantlly different phylumpro.", booktabs=TRUE)
#```
#
#(ref:taxonomy-phylumpro-$tmp[0]vs$tmp[1]) 热图展示$tmp[0]vs$tmp[1]在门及变形菌纲水平差异分类单元。Heatmap show differentially abundance phylumpro.[PDF](result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.pdf)
#
#```{r taxonomy-phylumpro-$tmp[0]vs$tmp[1], fig.cap="(ref:taxonomy-phylumpro-$tmp[0]vs$tmp[1])", out.width="99%"}
#knitr::include_graphics("result_k1-c/heat_phylumpro_$tmp[0]vs$tmp[1]_sig.png")
#```
#
#!;
#}else{
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1]
#
#无显著差异丰度分类单元；No significantlly differentially abundance taxonomy.
#
#!;
#}
#}
#
#
#
### 读group venn文件
#open DATABASE,"<$opts{v}";
#while (<DATABASE>) {
#	chomp;
#	push @venn,$_;
#}
#close DATABASE;
#
#if (@venn>=1) {
#
#print OUTPUT qq!
### 组间共有目Venn order {#result-order-venn}
#
#!;
#
#$i=0;
#foreach (@venn) {
#	chomp;
#	$i++;
#	my @tmp=split/\t/; # sampleA and sampleB
#	my $venn_list2;
#	foreach $tmp (@tmp) {
#		$venn_list2.=$tmp;
#	}
#	$tmp[2]="C" unless defined($tmp[2]);
#	$tmp[3]="D" unless defined($tmp[3]);
#	$tmp[4]="E" unless defined($tmp[4]);
#	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";
#
#print OUTPUT qq!
#
#### $venn_list
#
#(ref:order-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/order.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/order.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/order.txt.venn$venn_list2.xls.xls)  
#
#
#```{r order-venn-$i, fig.cap="(ref:order-venn-$i)", out.width="99%"}
#figs_2 = paste0("result_k1-c/order.txt.venn", "$venn_list", ".png")
#knitr::include_graphics(figs_2)
#```
#
#!;
#}
#
#print OUTPUT qq!
### 组间共有科Venn family  {#result-family-venn}
#
#!;
#$i=0;
#foreach (@venn) {
#	chomp;
#	$i++;
#	my @tmp=split/\t/; # sampleA and sampleB
#	my $venn_list2;
#	foreach $tmp (@tmp) {
#		$venn_list2.=$tmp;
#	}
#	$tmp[2]="C" unless defined($tmp[2]);
#	$tmp[3]="D" unless defined($tmp[3]);
#	$tmp[4]="E" unless defined($tmp[4]);
#	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";
#
#print OUTPUT qq!
#
#### $venn_list
#
#(ref:family-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/family.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/family.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/family.txt.venn$venn_list2.xls.xls)  
#
#
#```{r family-venn-$i, fig.cap="(ref:family-venn-$i)", out.width="99%"}
#figs_2 = paste0("result_k1-c/family.txt.venn", "$venn_list", ".png")
#knitr::include_graphics(figs_2)
#```
#
#!;
#}
#}
#
#print OUTPUT qq!
### 组间差异科分类学样式Pie family {#result-family-pie}
#
#(ref:family-pie) 比较组间差异科的门水平分类学样式。饼形图展示各种差异科细菌门水平分类比例。中间数字为所有显著差异科的数目，第一列为显著上调的科，第二列为显著下调的科，从上到下为各比较组。Pie charts show phylum of bacterial familys identified as either enriched or depleted in each genotype compared with WT. The number of familys in each category is noted inside each donut. !;
#
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#print OUTPUT qq!
#[$tmp[0]vs$tmp[1] enriched pie PDF](result_k1-c/pie_family_$tmp[0]vs$tmp[1]_enriched.pdf) 
#[$tmp[0]vs$tmp[1] depleted pie PDF](result_k1-c/pie_family_$tmp[0]vs$tmp[1]_depleted.pdf) !;
#$pie_list.="\"$tmp[0]vs$tmp[1]_enriched\"\, \"$tmp[0]vs$tmp[1]_depleted\"\, ";
#}
#$pie_list=~s/\,\ $//;
#
#print OUTPUT qq!
#
#```{r family-pie, fig.cap="(ref:family-pie)", out.width="49%"}
#figs_2 = paste0("result_k1-c/pie_family_", c(${pie_list}),".png")
#knitr::include_graphics(figs_2)
#```
#
#!;
#
#
#close OUTPUT;
#
#
#
#open OUTPUT,">$opts{o}04-f-otu.Rmd";
#print OUTPUT qq!
## OTUs差异分析 {#result-otu}
#
### 差异OTUs概述 {#result-otu-sum}
#
#样品组间显著差异OTUs数量(P < 0.05, FDR < 0.05)如 Table \\\@ref(tab:otu-sum) 所示。[TXT](result_k1-c/otu_sum.txt)；
#
#```{r otu-sum}
#table_otu = read.table("result_k1-c/otu_sum.txt", sep="\t", header=T)
#knitr::kable(table_otu, caption="各样品组间差异OTUs数量汇总", booktabs=TRUE)
#```
#
#
#样品组间显著差异显著差异OTU详细列表(P < 0.05, FDR < 0.05)如 Table \\\@ref(tab:otu) 所示。[TXT](result_k1-c/otu.txt)
#
#```{r otu}
#table_otu = read.table("result_k1-c/otu.txt", sep="\t", header=F)
#colnames(table_otu) = c("OTU","Sample A vs B","P-value")
#table_otu=head(table_otu,n=10)
#knitr::kable(table_otu, caption="样品组间显著差异前10个OTU, 完整表下载见上方TXT链接", booktabs=TRUE)
#```
#
#[下载OTU各组各样品相对丰度和物种信息 raw Data](result_k1-c/database.txt)
#
### 差异OTU {#result-otu-da}
#
#!;
#
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#	`ln -f result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.txt result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls`;
#	`ln -f result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.txt result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls`;
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1]
#
#(ref:otu-$tmp[0]vs$tmp[1]) $tmp[0]vs$tmp[1]基因型存在一些丰度显著上调或下调的OTU (P & FDR < 0.05, GLM likelihood rate test)。(A) 火山图展示$tmp[0]与$tmp[1]相比OTU的变化，x轴为OTU差异倍数取以2为底的对数，y轴为取丰度值百万比取2为底的对数，红蓝代表显著上下调；(B) 热图展示$tmp[0]与$tmp[1]显著差异OTU在每个样品中丰度值，数据采用Z-Score方法进行标准化，红色代表丰度相对高，而绿色代表丰度相对低，黄色代表中间水平；(C) 曼哈顿图展示OTU的变化情况及在各门水平中的分布，x轴为OTU按物种门水平物种注释字母排序，y轴为Pvalue值取自然对数，虚线为采用FDR校正的P-value的显著性阈值，图中每个点代表OTU，颜色为门水平注释，大小为相对丰度，形状为变化类型，其中上实心三角为显著上调，而下空心三角为显著下调；(D) 曼哈顿图按目水平上色。
#$tmp[0] are enriched and depleted for certain OTUs (P & FDR < 0.05, GLM likelihood rate test). (A) Volcano plot overview of abundance and fold change of OTUs; (B) Heatmap showing differentially abundance OTUs; (C) Manhattan plot showing phylum pattern of differentially abundance OTUs, colored by phylum; (D) Manhattan plot colored by order.
#[Volcano plot PDF](result_k1-c/vol_otu_$tmp[0]vs$tmp[1].pdf)  [Heatmap PDF](result_k1-c/heat_otu_$tmp[0]vs$tmp[1]_sig.pdf)  [Manhattan plot phlyum PDF](result_k1-c/man_otu_$tmp[0]vs$tmp[1].pdf) [Manhattan plot order PDF](result_k1-c/man_order_$tmp[0]vs$tmp[1].pdf)
#
#```{r otu-$tmp[0]vs$tmp[1], fig.cap="(ref:otu-$tmp[0]vs$tmp[1])", out.width="99%"}
#figs_2 = paste0("result_k1-c/", c("vol_otu_$tmp[0]vs$tmp[1]", "heat_otu_$tmp[0]vs$tmp[1]_sig", "man_otu_$tmp[0]vs$tmp[1]", "man_order_$tmp[0]vs$tmp[1]"),".png")
#knitr::include_graphics(figs_2)
#```
#
#$tmp[0]与$tmp[1]相比显著差异的OTU如 Table \\\@ref(tab:tab-otu-$tmp[0]vs$tmp[1]) 所示。[Enriched TXT](result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls)  [Depleted TXT](result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls)
#
#```{r tab-otu-$tmp[0]vs$tmp[1]}
#e = read.table("result_k1-c/otu_$tmp[0]vs$tmp[1]_enriched.xls", sep="\t", header=T)
#e = head(e[order(-e\$A_mean),],n=10)
#d = read.table("result_k1-c/otu_$tmp[0]vs$tmp[1]_depleted.xls", sep="\t", header=T)
#d = head(d[order(-d\$B_mean),],n=10)
#m = rbind(e,d)
#m=m[,c("otu","A_mean","B_mean","logFC","PValue","phylum","class","order","family","genus")]
#knitr::kable(m, caption="样品组$tmp[0] vs $tmp[1]相比显著差异的OTU详细信息；Significantlly different OTU.", booktabs=TRUE)
#```
#
#!;
#}
#
#print OUTPUT qq!
### 组间差异OTU分类学样式Pie OTU  {#result-otu-pie}
#
#(ref:otu-pie) 比较组间差异OTU的分类学样式。饼形图展示各种差异OTU细菌门水平分类比例，图例同上图曼哈顿图。中间数字为所有显著差异OTU的数目，第一列为显著上调的OTU，第二列为显著下调的OTU，从上到下为各比较组。Pie charts show phylum of bacterial OTUs identified as either enriched or depleted in each genotype compared with WT. The number of OTUs in each category is noted inside each donut. !;
#my $pie_list;
#my $venn_list;
#foreach (@group) {
#	chomp;
#	my @tmp=split/\t/; # sampleA and sampleB
#print OUTPUT qq!
#[$tmp[0]vs$tmp[1] enriched pie PDF](result_k1-c/pie_otu_$tmp[0]vs$tmp[1]_enriched.pdf) 
#[$tmp[0]vs$tmp[1] depleted pie PDF](result_k1-c/pie_otu_$tmp[0]vs$tmp[1]_depleted.pdf) !;
#$pie_list.="\"$tmp[0]vs$tmp[1]_enriched\"\, \"$tmp[0]vs$tmp[1]_depleted\"\, ";
#}
#$pie_list=~s/\,\ $//;
#
#print OUTPUT qq!
#
#```{r otu-pie, fig.cap="(ref:otu-pie)", out.width="49%"}
#figs_2 = paste0("result_k1-c/pie_otu_", c(${pie_list}),".png")
#knitr::include_graphics(figs_2)
#```
#
#!;
#
#if (@venn>=1) {
#print OUTPUT qq!
### 组间共有OTU(Venn)  {#result-otu-venn}
#
#!;
#
#my $i=0;
#foreach (@venn) {
#	chomp;
#	$i++;
#	my @tmp=split/\t/; # sampleA and sampleB
#	my $venn_list2;
#	foreach $tmp (@tmp) {
#		$venn_list2.=$tmp;
#	}
#	$tmp[2]="C" unless defined($tmp[2]);
#	$tmp[3]="D" unless defined($tmp[3]);
#	$tmp[4]="E" unless defined($tmp[4]);
#	$venn_list="$tmp[0]$tmp[1]$tmp[2]$tmp[3]$tmp[4]";
#
#print OUTPUT qq!
#
#### $venn_list
#
#(ref:otu-venn-$i) 维恩图展示各比较组差异OTU的共有和特有数量。Venn diagrams show common and unique OTUs in each group. [Figure PDF](result_k1-c/otu.txt.venn$venn_list.pdf)  [List XLS](result_k1-c/otu.txt.venn$venn_list2.xls)  [Detail XLSX](result_k1-c/otu.txt.venn$venn_list2.xls.xls)  
#
#
#```{r otu-venn-$i, fig.cap="(ref:otu-venn-$i)", out.width="99%"}
#figs_2 = paste0("result_k1-c/otu.txt.venn", "$venn_list", ".png")
#knitr::include_graphics(figs_2)
#```
#
#!;
#}
#
#}
#
##$file = "result_k1-c/ter_sum.txt";
#if (-e "doc/group_tern.txt") {
#print OUTPUT qq!
### 三元图Ternary 
#
#三元图中各部分高亮OTU数量统计信息如 Table \\\@ref(tab:ter-sum) 所示。[Summary table](result_k1-c/ter_sum.txt)
#
#```{r ter-sum}
#tab = read.table("result_k1-c/ter_sum.txt", sep="\\t", header=T)
#knitr::kable(tab, caption="三元图中各部分高亮OTU数量统计", booktabs=TRUE)
#```
#!;
#
#open DATABASE,"<$opts{t}";
#while (<DATABASE>) {
#	chomp;
#	@tmp=split/\t/;;
#	foreach (@tmp) {
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.txt > result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.xls` if -e "result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$_.txt";
#	`sed ':a;N;\$!ba;s/^/ID\\t/g' result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.txt > result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.xls` if -e "result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]${_}venn.txt";
#	}
#print OUTPUT qq!
#### $tmp[0] vs $tmp[1] vs $tmp[2]
#
#(ref:ter-$tmp[0]$tmp[1]$tmp[2]) 三元图展示$tmp[0]、$tmp[1]、和$tmp[2]三组OTU(>1‰)的相对丰度。(A) 各组相对于另外两组显著差异的OTU  [PDF](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2].pdf)。每个点代表一个OTU，位置代表它在各组间的相对比例，大小代表三组的平均丰度。彩色的圆圈代表每组中相对于其它两组显著富集的OTU，位于三角形[左$tmp[0]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[0].xls)、[右$tmp[1]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[1].xls)和[顶$tmp[2]](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[2].xls)组所特异显著富集的OTU，分别用绿、橙和红高亮显示；
#(B) 三角形左$tmp[0]、右$tmp[1]相对于顶$tmp[2] 组富集的OTU  [PDF](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]venn.pdf)。图中底部[两组共有显著富集OTU用红色显示](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[2]venn.xls)，[左$tmp[0]特异的OTU分别用绿色显示](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[0]venn.xls)，[右组特异显示为橙色](result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]$tmp[1].xls)。
#Ternary plot depicting group relative abundances (RAs) of all OTUs (>1‰) for $tmp[0], $tmp[1], and $tmp[2]. (A) Compared with each others, three groups specific OTU. Each point corresponds to an OTU. Its position represents its RA with respect to each group, and its size represents the average across all three groups. Colored circles represent OTUs enriched in one group compared with the others (green in $tmp[0], orange in $tmp[1], and red in $tmp[2] samples), whereas gray circles represent OTUs that are not significantly enriched in a specific groups. (B) Group $tmp[2] (top) as control, $tmp[1] and $tmp[2] specific and common OTU; Colored circles represent OTUs enriched in $tmp[0] (green), $tmp[1] (orange), and common in $tmp[0] and $tmp[1] (red).
#
#```{r ter-$tmp[0]$tmp[1]$tmp[2], fig.cap="(ref:ter-$tmp[0]$tmp[1]$tmp[2])", out.width="99%"}
#figs_1 = paste0("result_k1-c/ter_$tmp[0]$tmp[1]$tmp[2]", c("", "venn"),".png")
#knitr::include_graphics(figs_1)
#```
#
#!;
#}
#}
#close OUTPUT;



open OUTPUT,">$opts{o}99-references.Rmd";
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
`rm -f $opts{b}/$opts{b}`;
`ln /mnt/bai/yongxin/test/usearchPE250/$opts{b} /var/www/html/report/16s/$opts{b} -sf`;

print "Result please visiting http://bailab.genetics.ac.cn/report/16s/$opts{b}\n";
