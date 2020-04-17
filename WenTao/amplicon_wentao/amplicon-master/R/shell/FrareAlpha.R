#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：Alpha稀释曲线绘制
# Functions: Alpha rarefaction curve

options(warn = -1) # Turn off warning


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
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小
#
# 抽平最小序列数-b，“begin”：最少抽平100条序列。
#
# 间隔抽平序列 -s，“step”：设置用于间隔多少条徐序列进行抽平，默认100条
#
#稀释曲线算法 -m ：可选用一下多种算法：observed , chao1  , diversity_inverse_simpson , diversity_gini_simpson,
# diversity_shannon   ,   diversity_fisher   ,  diversity_coverage     ,    evenness_camargo,
# evenness_pielou    ,   evenness_simpson       ,    evenness_evar ,   evenness_bulla,
# dominance_dbp      ,  dominance_dmn        ,      dominance_absolute   ,      dominance_relative,
# dominance_simpson      ,    dominance_core_abundance ,  dominance_gini  ,           rarity_log_modulo_skewness,
# rarity_low_abundance   ,    rarity_noncore_abundance,  rarity_rare_abundance






site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="./otutab.txt",
                help="Alpha rarefaction richness [default %default]"),
    make_option(c("-d", "--design"), type="character", default="./metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=170,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=110,
                help="Figure heidth in mm [default %default]"),
    make_option(c("-b", "--begin"), type="numeric", default=100,
                help="rera begin sample number [default %default]"),
    make_option(c("-s", "--step"), type="numeric", default=100,
                help="rera interval [default %default]"),
    make_option(c("-m", "--method"), type="character", default="chao1",
                help="method used calculate reraplot:observed , chao1  , diversity_inverse_simpson , diversity_gini_simpson,
        diversity_shannon   ,   diversity_fisher   ,  diversity_coverage     ,    evenness_camargo,
        evenness_pielou    ,   evenness_simpson       ,    evenness_evar ,   evenness_bulla,
        dominance_dbp      ,  dominance_dmn        ,      dominance_absolute   ,      dominance_relative,
        dominance_simpson      ,    dominance_core_abundance ,  dominance_gini  ,           rarity_log_modulo_skewness,
        rarity_low_abundance   ,    rarity_noncore_abundance,  rarity_rare_abundance [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}


# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output="./"}




#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
map = read.delim(opts$design,row.names = 1)
#----2.2 otu table otu#----
otu = read.delim(opts$input,row.names = 1)


source("../rareAlphaplot-wt.R")

result1 = rareAlpha(otu = otu,map = map,group = opts$group, start =opts$begin,step = opts$step,method = opts$method)
# 提取图形
p1 =  result1[[1]]

# 提取数据
data = result1[[2]]
# head(data)
# 提取分组稀释曲线
p2 =  result1[[3]]

# 提取分组稀释曲线
p3 =  result1[[4]]



#---3.2 保存 Saving#----
# 大家可以修改图片名称和位置，长宽单位为毫米
library(ggplot2)
o=paste0(opts$output,"/rareplot.pdf")
ggsave(o, p1, width = opts$width, height = opts$height, units = "mm")

o=paste0(opts$output,"/rareGroupPlot.pdf")
ggsave(o, p2, width = opts$width, height = opts$height, units = "mm")
o=paste0(opts$output,"/rareGroupBarPlot.pdf")
ggsave(o, p3, width = opts$width, height = opts$height, units = "mm")
o=paste0(opts$output,"/rareplotData.csv")
write.csv(data,o)


# ggsave(paste0(opts$output,"/alpha_boxplot_",opts$alpha_index,".pdf"), p, width = opts$width, height = opts$height, units = "mm")




