# 绘制韦恩图和upset图表
#
# This is the first function named 'VenUpset'
# which draw Ven plot and Upset plot with otutab and metadata, and reture  a base plot object
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

#' @title Plotting Ven and Upset plot for each group
#' @description Input otutab and metadata
#' VennDiagram::Ven plot
#' UpSetR::Upset plot
#' @param otutab rarefied OTU table, typical output of usearch -otutab_norm or -otutab_rare,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param group column name for groupID.
#' @param rep repeat number of each group
#' @return base plot object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso Ven-Upset
#' @examples
#' # Set four parameters: otutab, metadata, group and rep
#' VenUpset(otu = otutab,map = metadata,group = "genotype",rep = 6)
#' @export

#代码测试


# #清空内存
# rm(list=ls())
# load("../../data/otutab.rda")
# load("../../data/metadata.rda")

#
## upset 我不会直接输出一个图形对象，所以这里直接运行就出结果，不知道行不行
#ven图可以提取出来，使用grid.draw出图。
# result = VenUpset(otu = otutab,map = metadata,group = "genotype",rep = 6)
# ven 出图
# grid.draw(result[[1]])
# 
# 
# otu = otutab
# map = metadata
# group = "genotype"
# rep = 6
# path = ""

VenUpset = function(otu = otutab,map = metadata,group = "genotype",rep = 6,path = ""){
  path = path
  # path = "./VenUpset/"
  # dir.create(path)
  library("phyloseq")
  library(UpSetR)
  # library("tibble")
  library (VennDiagram) 
  
  # 依赖关系检测与安装
  p_list = c("phyloseq", "dplyr", "UpSetR","VennDiagram") 
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  }
  
  
  #导入otu表格
  otu = otu
  head(otu)
  otu = as.matrix(otu)


  #导入分组文件
  map = metadata
  head(map)

  
  ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                 sample_data(map) 
                 # phy_tree(tree)
                 
  )
  
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  aa = vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  # dim(count)
  sub_design <- as.data.frame(sample_data(ps))
  
  # levels(sub_design$SampleType)[1]
  # name1 = paste("name",levels(sub_design$SampleType)[1],sep = "")
  # sub_design $SampleType
  # levels(sub_design $SampleType)
  ##########这里的操作为提取三个分组
  pick_val_num <- rep/2
  count[count > 0] <- 1###这个函数只能用于0,1 的数据，所以我这么转换
  
  count2 = as.data.frame(count)
  
  
  aa = sub_design[,group]
  colnames(aa) = "Vengroup"
  
  #数据分组
  iris.split <- split(count2,as.factor(aa$Vengroup))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  # head(ven2)
  ven2[ven2 < pick_val_num]  = 0
  ven2[ven2 >=pick_val_num]  = 1
  ven2 = as.data.frame(ven2)
  
  
  #########更加高级的设置在这里可以查看#https://mp.weixin.qq.com/s/6l7gftKQfiyxNH66i19YtA
  ven3 = as.list(ven2)
  
  ven_pick = get.venn.partitions(ven3)
  
  for (i in 1:ncol(ven2)) {
    
    
    ven3[[i]] <-  row.names(ven2[ven2[i] == 1,])
    
  }
  
  
  if (length(names(ven3)) == 2) {
    filename3 = paste(path,"/ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename3,width = 8, height = 6)
    T<-venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue"), #填充颜色
                    col=c('red',"blue"), #圈线颜色
                    cat.col=c('red',"blue"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    # grid.draw(T)
    # dev.off();
    
    
    
    
    
  } else if (length(names(ven3)) == 3) {
    filename3 = paste(path,"/ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename3,width = 12, height = 12)
    T<-venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow"), #填充颜色
                    col=c('red',"blue","yellow"), #圈线颜色
                    cat.col=c('red',"blue","yellow"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    # grid.draw(T)
    # dev.off()
    grid.draw(T)
  } else if (length(names(ven3)) == 4) {
    filename3 = paste(path,"/ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename3,width = 12, height = 12)
    T<-venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow","#7ad2f6"), #填充颜色
                    col=c('red',"blue","yellow","#7ad2f6"), #圈线颜色
                    cat.col=c('red',"blue","yellow","#7ad2f6"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    # grid.draw(T)
    # dev.off()
    grid.draw(T)
  }else if (length(names(ven3)) == 5) {
    filename3 = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename3,width = 12, height = 12)
    T<-venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow","#7ad2f6","green"), #填充颜色
                    col=c('red',"blue","yellow","#7ad2f6","green"), #圈线颜色
                    cat.col=c('red',"blue","yellow","#7ad2f6","green"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    # grid.draw(T)
    # dev.off()
    grid.draw(T)
  }else if (length(names(ven3)) == 6) {
    
    print("ven not use for more than 6")
  }
  
  dev.off()
  
  
  if (ven_pick[[1,6]] != 0) {
    
    filename4 = paste(path,"./UpSet_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename4,width = 12, height = 8)
    upset(ven2, sets = colnames(ven2),
          number.angles = 30, point.size = 2, line.size = 1,
          mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
          text.scale = c(2, 2, 2,2, 2, 2),mb.ratio = c(0.7, 0.3),order.by = "freq",keep.order = TRUE,
          queries = list(list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T)))
    
    
    dev.off()
    
    upset(ven2, sets = colnames(ven2),
          number.angles = 30, point.size = 2, line.size = 1,
          mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
          text.scale = c(2, 2, 2,2, 2, 2),mb.ratio = c(0.7, 0.3),order.by = "freq",keep.order = TRUE,
          queries = list(list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T),
                         list(query = intersects, params = 
                                list(colnames(ven2)), color = "red", active = T)))
  }
 return(list(T,ven2))
}


