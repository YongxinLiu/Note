
#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：对微生物不同门类挑选丰度较大的进行展示
# Functions: beta diversity plot

## 参数解释
## 设置输入输出文件和参数
# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始otu表格(otutab.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，result/alpha/vegan.txt; alpha多样性表格
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-o", "--output"，默认为输出目录，图片文件名为alpha_boxplot_+多样性指数名+.pdf；统计文本位于代码运行目录中alpha_boxplot_TukeyHSD.txt；
#
# 图片宽"-w", "--width"，默认200 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认140 mm，根据图像布局可适当增大或缩小
#
#"-N", "--maxnum"，默认10，选择丰度最高的前十个物种进行展示；
#
#"-T", "--TaxRank"，默认Phylum，选择需要展示的物种门类水平


options(warn = -1) # Turn off warning


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="otutab.txt",
                help="otu table [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-t", "--tax"), type="character", default="taxonomy.txt",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=200,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=140,
                help="Figure heidth in mm [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-N", "--maxnum"), type="numeric", default=10,
                help="Group name [default %default]"),
    make_option(c("-T", "--TaxRank"), type="character", default="Phylum",
                help="Group name [default %default]")


  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}

source("../BubblePlot-wt.R")

# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf

if(opts$output==""){opts$output="./"}

# opts$design = "../data/metadata.tsv"
# 
# opts$input = "../data/otutab.txt"
# 
# opts$tax ="../data/taxonomy.txt"
#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
map = read.delim(opts$design,row.names = 1)
#----2.2 otu table otu#----
otu = read.delim(opts$input,row.names = 1)

# tree = read_tree("./otus.tree")



tax = read.delim(opts$tax,row.names = 1)



result = BubblePlot(otu = otu,tax = tax,map = map ,N = opts$maxnum,taxGlomRank = opts$TaxRank,group = opts$group)

# 提取气泡图:样本为纵向
p1 = result[[1]]
# 样本为横向
p2 = result[[2]]
# p2


#提取作图数据
data1 = result[[3]]


o=paste0(opts$output,paste("/","BubbleplotR.pdf",sep = ""))
ggsave(o, p1, width = opts$width, height = opts$height, units = "mm")
o=paste0(opts$output,paste("/","BubbleplotC.pdf",sep = ""))
ggsave(o, p2, width = opts$height, height = opts$width, units = "mm")


o=paste0(opts$output,"/Bubbleplot.csv")
write.csv(data1,o)
