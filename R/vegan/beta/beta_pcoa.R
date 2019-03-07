#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1. 程序功能描述和主要步骤

# 程序功能：Beta多样性主坐标轴分析及组间统计
# Functions: PCoA analysis of samples and groups comparing
# Main steps: 
# - Reads distance matrix input.txt
# - Calculate orrdinate by PCoA and show in scatter plot
# - Adonis calculate significant between groups distance and group inner distance

# 程序使用示例
# USAGE
# # 展示样品间距离分布，统计组间是否显著，也用于异常样品筛选
# 
# Rscript ./script/beta_pcoa.r -h
# 
# # 默认基于bray_curtis距离
# Rscript ./script/beta_pcoa.r
# 
# # 完整默认参数
# Rscript ./script/beta_pcoa.r -i beta/bray_curtis.txt -t bray_curtis \
# -d doc/design.txt -n group \
# -o beta/pcoa_bray_curtis \
# -w 4 -e 2.5 
# 
# # 基于unifrac距离
# Rscript ./script/beta_pcoa.r -t unifrac
options(warn = -1)



# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}
# 解析命令行
if (TRUE){
  option_list = list(
    make_option(c("-t", "--type"), type="character", default="bray_curtis",
                help="Distance type; 距离类型, 可选bray_curtis, bray_curtis_binary, euclidean, jaccard, jaccard_binary, manhatten, unifrac, unifrac_binary [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="",
                help="Input beta distance; 距离矩阵,默认beta目录下与t同名，可指定 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="../design.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有txt和矢量图pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("",opts$type, ".txt", sep = "")}
  if (opts$output==""){opts$output=paste("pcoa_",opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The distance matrix file is ", opts$input,  sep = ""))
  print(paste("Type of distance type is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("digest","ggrepel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list = c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}




# 3. 读取输入文件

# 读取距离矩阵文件
dis = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 

# 提取样品组信息,默认为group可指定
sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
colnames(sampFile)[1] = "group"



# 4. 统计与绘图

# vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前3维
pcoa = cmdscale(dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
eig = pcoa$eig
points = cbind(points, sampFile[rownames(points),])
colnames(points) = c("x", "y", "z","group") 

# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group)) + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=paste(opts$type," PCoA",sep=""))  + 
  stat_ellipse(level=0.68) + theme_classic()
p
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, ".pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, ".pdf finished.", sep = ""))



# 添加样品标签
p=p+geom_text_repel(label=paste(rownames(points)),colour="black",size=3)
p
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_label.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_label.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_label.pdf finished.", sep = ""))


# Compare each group beta by vegan adonis in bray_curtis
da_adonis = function(sampleV){
  sampleA = as.matrix(sampleV$sampA)
  sampleB = as.matrix(sampleV$sampB)
  design2 = subset(sampFile, group %in% c(sampleA,sampleB))
  if (length(unique(design2$group))>1) {
    sub_dis_table = dis_table[rownames(design2),rownames(design2)]
    sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    print(paste("In ",opts$type," pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
    adonis_pvalue = paste(opts$type, sampleA, sampleB, adonis_pvalue, sep="\t")
    write.table(adonis_pvalue, file=paste(opts$output, ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
  }
}

# loop for each group pair
dis_table = as.matrix(dis)
if (TRUE) {
  compare_data = as.vector(unique(design[[opts$group]]))
  len_compare_data = length(compare_data)
  for(i in 1:(len_compare_data-1)) {
    for(j in (i+1):len_compare_data) {
      tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
      print(tmp_compare)
      da_adonis(tmp_compare)
    }
  }
}else {
  compare_data = read.table("doc/compare.txt", sep="\t", check.names=F, quote='', comment.char="")
  colnames(compare_data) = c("sampA", "sampB")
  for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	 
print(paste("Adnois statistics result in",opts$output, ".txt is finished.", sep = ""))

# 5. 保存图表

# 提示工作完成
print(paste("Output in ", opts$output, ".txt/pdf finished.", sep = ""))

