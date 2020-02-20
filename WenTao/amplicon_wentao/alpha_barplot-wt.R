# 绘制alpha多样性箱线图并添加统计分组 Alpha boxplot + LSD.test
#
# This is the first function named 'alpha_barplot'
# which draw boxplot with alpha and metadata, and reture a ggplot2 object
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate doc:              'Ctrl + Shift + Alt + R'

#' @title Plotting alpha diversity boxplot for each group with anova statistics
#' @description Input alpha index and metadata, and manual set alpha index and metadata column names.
#' agricolae::LSD.test calculate p-value, and dplyr summary each group max for p-value groups position.
#' ggplot2 show boxplot, jitter and stat groups.
#' @param alpha_div alpha diversity matrix, typical output of usearch -alpha_div,
#' rowname is sampleID, colname is index of alpha diversity;
#' @param index index(type) of alpha diversity;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @details
#' By default, returns richness diversity index
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: chao1, richness, shannon_e}
#' \item{other used indices: berger_parker, buzas_gibson, dominance, equitability, jost, jost1, reads, robbins, simpson, shannon_2, shannon_10}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso alpha_rare
#' @examples
#' # Set four parameters: alpha_div, metadata, index and groupID
#' alpha_barplot(alpha_div, metadata, "richness", "genotype")
#' # Set two parameters: alpha_div, metadata, and index and groupID as default richness and genotype
#' alpha_barplot(alpha_div, metadata)
#' # Set two parameters: alpha_div, metadata, and index and groupID as using chao1 and site
#' alpha_barplot(alpha_div, metadata, "chao1", "site")
#' @export

#代码测试
# load("../data/alpha_div.rda")
# load("../data/metadata.rda")
# p = alpha_barplot(alpha_div, metadata, index = "richness", groupID = "genotype")
# p
# index = "richness"
# groupID = "genotype"



alpha_barplot <- function(alpha_div, metadata, index = "richness", groupID = "genotype") {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "dplyr", "multcompView") # "agricolae"
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  }

  # 测试默认参数
  # library(amplicon)
  # index = "richness"
  # groupID = "genotype"
  # metadata = subset(metadata, genotype %in% c("KO","OE"))

  # 交叉筛选
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx,]
  alpha_div = alpha_div[rownames(metadata),]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  # colnames(sampFile)[1] = "group"

  # 合并alpha_div和metadata
  df = cbind(alpha_div[rownames(sampFile),index], sampFile)
  colnames(df) = c(index,"group")
  df
  
  
  # 统计各种显著性
  model = aov(df[[index]] ~ group, data=df)
  # 计算Tukey显著性差异检验
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  # 提取比较结果
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)

  # 保存统计结果
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table(paste(date(), "\nGroup\t", groupID, "\n\t", sep=""), file=paste("alpha_boxplot_TukeyHSD.txt",sep=""),append = T, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  suppressWarnings(write.table(Tukey_HSD_table, file=paste("alpha_boxplot_TukeyHSD.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

  # 函数：将Tukey检验结果P值转换为显著字母分组
  # 输入文件为图基检验结果和分组
  generate_label_df = function(TUKEY, variable){
    # library(multcompView)
    # 转换P值为字母分组
    ## 提取图基检验中分组子表的第4列P adjust值
    Tukey.levels = TUKEY[[variable]][,4]
    ## multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
    Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])

    # 按分组名字母顺序
    ## 提取字母分组行名为group组名
    Tukey.labels$group = rownames(Tukey.labels)
    # 按组名的字母顺序排列，默认的Levels
    Tukey.labels=Tukey.labels[order(Tukey.labels$group), ]
    return(Tukey.labels)
  }

  # 当只有两组时，用LSD标注字母
  if (length(unique(df$group)) == 2){
    # LSD检验，添加差异组字母
    library(agricolae)
    out = LSD.test(model, "group", p.adj="none")
    stat = out$groups
    # 分组结果添入Index
    df$stat=stat[as.character(df$group),]$groups
  # 当大于两组时，用multcompView标注字母
  }else{
    # library(multcompView)
    LABELS = generate_label_df(Tukey_HSD , "group")
    df$stat=LABELS[as.character(df$group),]$Letters
  }

  # 设置分组位置为各组y最大值+高的5%
  max=max(df[,c(index)])
  min=min(df[,index])
  x = df[,c("group",index)]
  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',index,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$group
  df$y=y[as.character(df$group),]$Max + (max-min)*0.05

  head(df)
  #求取均值和方差
  wen1 = as.data.frame(tapply(as.vector(as.matrix(df[1])),df$group,mean,na.rm=TRUE))
  wen2 = as.data.frame(tapply(as.vector(as.matrix(df[1])),df$group,sd,na.rm=TRUE))
  went = cbind(wen1,wen2)
  
  colnames(went) = c("mean" ,"SD")
  went
  #去除重复，用来提取显著性字母
  aa = distinct(df, group, .keep_all = TRUE)  
  #将显著性字母添加到表格中
  went$label = aa$stat[match(row.names(went),aa$group)] 
  went$ID = row.names(went)
  #显著性标记由于坐标轴长度不够会遮盖字母，所以我这里人工设定一个y轴坐标大小，是最大值的1.1倍
  a = max(went$mean + went$SD)*1.1
  
  # 绘图 plotting

  p = ggplot(went , aes(x = ID, y = mean,colour= ID)) +
    geom_bar(aes(colour= ID,fill = ID),stat = "identity", width = 0.4,position = "dodge") +
    scale_y_continuous(expand = c(0,0),limits = c(0,a)) +
    geom_errorbar(aes(ymin=mean-SD,
                      ymax=mean+SD),
                  colour="black",width=0.1,size = 1)+
    geom_text(aes(label = label,y=mean+SD, x = ID,vjust = -0.3),color = "black")+theme_classic() +
    theme(text=element_text(family="sans", size=7))
  p
  
}
























