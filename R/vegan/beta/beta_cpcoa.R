#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>



# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory 设置工作目录为 data (分析项目根目录)
# 清空工作环境 Clean enviroment object
rm(list=ls()) 



# 1. 程序功能描述和主要步骤

# 程序功能：限制性主坐标轴分析及组间统计
# Functions: Constrained PCoA analysis of groups
# Main steps: 
# - Reads OTU table result/otutab.txt
# - Orrdinate by CCA and show in scatter plot
# - Aov.cca() calculate significant between all groups distance

# 程序使用示例
# # 显示脚本帮助
# Rscript ../6script/beta_cpcoa.R -h
# 
# # 基于bray距离计算CCA，默认bray方法
# Rscript ../6script/beta_cpcoa.R
# 
# # 基于jaccard距离计算CCA
# Rscript ../6script/beta_cpcoa.R -t jaccard
# 
# # 附完整参数
# Rscript ../6script/beta_cpcoa.R -i gg/otutab.txt -t bray \
# -d doc/design.txt -n group \
# -o gg/beta/cpcoa_bray # filename prefix for output directory name 
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
    make_option(c("-t", "--type"), type="character", default="bray",
                help="Distance type; 距离类型, 可选manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup , binomial, chao, cao or mahalanobis[default %default]"),   
    make_option(c("-i", "--input"), type="character", default="otutab_norm.txt",
                help="Input norm otutab; 标准化OTU表 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="../design.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 通常会有统计表txt、矢量图pdf和位图png [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){opts$output=paste("cpcoa_",opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The normalized OTU table file is ", opts$input,  sep = ""))
  print(paste("Distance type is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}



# 2. 依赖关系检查、安装和加载

# 2.1 安装CRAN来源常用包
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
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 

# 提取样品组信息,默认为group可指定
sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
colnames(sampFile)[1] = "group"
sampFile$sample=row.names(design)

# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(otutab) # match design with alpha
sampFile = sampFile[idx,]
otutab = otutab[,rownames(sampFile)] 



# 4. 统计与绘图

# 提取CCA中主要结果
variability_table = function(cca){
  chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table = cbind(chi, chi/chi[1])
  colnames(variability_table) = c("inertia", "proportion")
  rownames(variability_table) = c("total", "constrained", "unconstrained")
  return(variability_table)
}

# Constrained analysis OTU table by genotype
capscale.gen = capscale(t(otutab) ~ group, data=sampFile, add=F, sqrt.dist=T, distance=opts$type) 

# ANOVA-like permutation analysis
perm_anova.gen = anova.cca(capscale.gen, permutations = 1000, parallel = 4)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen = variability_table(capscale.gen)
eig = capscale.gen$CCA$eig
variance = var_tbl.gen["constrained", "proportion"]
p.val = perm_anova.gen[1, 4]

# extract the weighted average (sample) scores
points = capscale.gen$CCA$wa[, 1:2]
points = as.data.frame(points)
colnames(points) = c("x", "y")
points = cbind(points, sampFile[match(rownames(points), rownames(sampFile)),])

# plot CPCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=group)) + geom_point() +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) + 
  theme_classic()+ 
  stat_ellipse(level=0.68)
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


# 5. 保存图表

# 提示工作完成

print("Beta diversity: Constrianed PCoA / CCA done!!!")
