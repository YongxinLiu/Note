
#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：beta多样性排序 三种方法群落差异检测，群落样本分组两组之间检测
# Functions: beta diversity plot

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
    # make_option(c("-i", "--input"), type="character", default="../otutab.txt",
    #             help="Alpha rarefaction richness [default %default]"),
    # make_option(c("-d", "--design"), type="character", default="../metadata.tsv",
    #             help="Design file or metadata [default %default]"),
    make_option(c("-i", "--input"), type="character", default="ps_liu.rds",
                help="Alpha rarefaction richness [default %default]"),
    # make_option(c("-xx", "--input"), type="character", default="ps_liu.rds",
    #             help="Community difference detection method [default %default]"),
    
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    # make_option(c("-w", "--width"), type="numeric", default=170,
    #             help="Figure width in mm [default %default]"),
    # make_option(c("-e", "--height"), type="numeric", default=110,
    #             help="Figure heidth in mm [default %default]"),
    # make_option(c("-m", "--method"), type="character", default="DCA",
    #             help="method could chioce: DCA, CCA, RDA, NMDS, MDS, PCoA,PCA,LDA,t-sne [default %default]"),
    make_option(c("-p", "--permutation"), type="character", default="adonis",
                help="Community difference detection method [default %default]"),
    make_option(c("-x", "--distance"), type="numeric", default=8,
                help="Community difference detection method [default %default]"),
    # make_option(c("-t", "--tree"), type="character", default="",
    #             help="Community difference detection method [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]")
    
    
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}
# opts$design
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf

# if(opts$output==""){opts$output="./"}
opts$input

#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
# map = read.delim(opts$design,row.names = 1)
# #----2.2 otu table otu#----
# otu = read.delim(opts$input,row.names = 1)
#   otu = phyloseq::otu_table(otu,taxa_are_rows=TRUE)
#   map = phyloseq::sample_data(map)
#   ps = phyloseq::merge_phyloseq(otu,map)
ps = readRDS(opts$input)

# tree = read_tree("./otus.tree")

# source("./BetaDiv-wt.R")
# source("./pairMicroTest-wt.R")
source("./pairMicroTest-wt.R")

# ps = readRDS("./ps_liu.rds")
# opts$distance
dist = opts$distance
# print(dist)
dist_methods <- unlist(phyloseq::distanceMethodList)
# print(dist_methods)
# print(dist)
dist1 =   dist_methods[dist]
# print(dist1)
# unif <- phyloseq::distance(ps, method=dist1)


result =pairMicroTest (ps = ps,dist = opts$distance,Micromet = opts$permutation)


#提取结果

o=paste0(opts$output,"/betadiff.csv")
write.csv(result,o)



