# library("vegan")
# library("grid")
# library("gridExtra")
# library(phyloseq)


# 全部的排序方法
# # methodAll = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA","PCA","LDA","t-sne")

# otu = read.delim("../otutab.txt",row.names = 1)
# map = read.delim("../metadata.tsv",row.names = 1)
# tree = read_tree("./otus.tree")
# result = BetaDiv(otu,map,dist = "bray",method ="DCA",group = "Group",pvalue.cutoff = 0.05,Micromet = "MRPP")
# result = BetaDiv(otu,map,tree,dist = "unifrac",method ="NMDS",pvalue.cutoff = 0.05,Micromet ="anosim")
# result = BetaDiv(otu,map,dist = "bray",method ="t-sne",pvalue.cutoff = 0.05,Micromet ="adonis")
# 
# ps = readRDS("./ps_liu.rds")
# result = BetaDiv(ps = ps,dist = "bray",method ="DCA",pvalue.cutoff = 0.05,Micromet = "MRPP")
# result = BetaDiv(ps = ps,dist = "unifrac",method ="NMDS",pvalue.cutoff = 0.05,Micromet ="anosim")
# result = BetaDiv(ps = ps,dist = "bray",method ="t-sne",pvalue.cutoff = 0.05,Micromet ="adonis")
# otu
# map
# tree = tree
# dist = opts$distance
# group = opts$group
# method =opts$method
# pvalue.cutoff = 0.05
# Micromet = opts$permutation
# 


BetaDiv = function(otu = NULL,map = NULL,tree = NULL,ps = NULL,group = "Group",dist = "bray",method ="DCA",
                   Micromet = "MCPP",pvalue.cutoff = 0.05){
  
  dist_methods <- unlist(phyloseq::distanceMethodList)
  dist =   dist_methods[dist]
  dist
  # 需要的R包
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  
  # 数据导入
  if (is.null(otu)&is.null(map)) {
    ps = ps
  }else {
    #导入otu表格
    otu = as.matrix(otu)
    str(otu)
    # #导入注释文件
    # tax = as.matrix(tax)
    # taxa_names(tax)
    
    #导入分组文件
    colnames(map) = gsub(group,"AA", colnames(map))
    
    map$Group = map$AA
    map$Group = as.factor(map$Group )
    map$Group 
    # #导入进化树
    # tree = read.tree("./otus.tree")
    # tree
    ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                   sample_data(map) 
                   # phy_tree(tree)
    )
  }
  
  # dist_methods <- unlist(distanceMethodList)
  if (dist %in% c("unifrac" , "wunifrac",  "dpcoa")) {
    phy_tree(ps) = tree
  }
  


# 求取相对丰度
ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela 

# 先择方法进行运算
if (method == "DCA") {
  #---------如果选用DCA排序
  ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
  #样本坐标
  points = ordi$rproj[,1:2]
  colnames(points) = c("x", "y") #命名行名
  #提取特征值
  eig = ordi$evals^2
}


# method = "CCA"
#---------如果选用CCA排序
if (method == "CCA") {
# dist = "bray"
ordi = ordinate(ps1_rela, method=method, distance=dist)
#样本坐标,这里可选u或者v矩阵
points = ordi$CA$v[,1:2]
colnames(points) = c("x", "y") #命名行名
#提取特征值
eig = ordi$CA$eig^2
}
# method ="RDA"
#---------如果选用RDA排序
if (method == "RDA") {
  # dist = "bray"
  ordi = ordinate(ps1_rela, method=method, distance=dist)
  #样本坐标,这里可选u或者v矩阵
  points = ordi$CA$v[,1:2]
  colnames(points) = c("x", "y") #命名行名
  #提取特征值
  eig = ordi$CA$eig
}

#---------如果选用DPCoA排序 不用做了，不选择这种方法了，这种方法运行太慢了
# method = "MDS"
#---------如果选用MDS排序 
if (method == "MDS") {
  # dist = "bray"
  ordi = ordinate(ps1_rela, method=ord_meths[i], distance=dist)
  #样本坐标,
  points = ordi$vectors[,1:2]
  colnames(points) = c("x", "y") #命名行名
  #提取解释度
  eig = ordi$values[,1]
}

# method = "PCoA"
#---------如果选用pcoa排序 
if (method == "PCoA") {
  # dist = "bray"
  unif <- phyloseq::distance(ps1_rela , method=dist, type="samples")
  #这里请记住pcoa函数
  pcoa = cmdscale(unif, k=2, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  
  points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
  colnames(points) = c("x", "y") #命名行名
  eig = pcoa$eig
}

#---------------------------------------PCA分析
# method = "PCA"
vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}

otu_table = as.data.frame(t(vegan_otu(ps1_rela )))
# head(otu_table)
if (method == "PCA") {
  
  otu.pca <- prcomp(t(otu_table), scale. = TRUE)
  #提取坐标
  points = otu.pca$x[,1:2]
  colnames(points) = c("x", "y") #命名行名
  # #提取荷载坐标
  # otu.pca$rotation
  # 提取解释度,这里提供的并不是特征值而是标准差，需要求其平方才是特征值
  eig=otu.pca$sdev
  eig=eig*eig
}
# method = "LDA"

#---------------LDA-
if (method == "LDA") {
  #拟合模型
  library(MASS)
  data = t(otu_table)
  # head(data)
  data = as.data.frame(data)
  # data$ID = row.names(data)
  # 
  data <- scale(data, center = TRUE, scale = TRUE)
  
  dim(data)
  data1 = data[,1:10]
  map = as.data.frame(sample_data(ps1_rela))
  
  model <- lda(data, map$Group)

  # 提取坐标
  ord_in = model
  axes = c(1:2)
  points <- data.frame(predict(ord_in)$x[, axes])
  colnames(points) = c("x", "y") #命名行名
  # 提取解释度
  eig= ord_in$svd^2

}

# method = "NMDS"
if (method == "NMDS") {
  #---------如果选用NMDS排序 
  # i = 5
  # dist = "bray"
  ordi = ordinate(ps1_rela, method=method, distance=dist)
  #样本坐标,
  points <- ordi$points[,1:2]
  colnames(points) = c("x", "y") #命名行名
  #提取stress
  stress = ordi$stress
  stress= paste("stress",":",round(stress,2),sep = "")
}



if (method == "t-sne") {
  data = t(otu_table)
  # head(data)
  data = as.data.frame(data)
  # data$ID = row.names(data)
  # 
  data <- scale(data, center = TRUE, scale = TRUE)
  
  dim(data)
  map = as.data.frame(sample_data(ps1_rela))
  row.names(map)
  #---------tsne
  # install.packages("Rtsne")
  library(Rtsne)
  
  tsne <- Rtsne(data,perplexity = 3)
  
  # 提取坐标
  points = as.data.frame(tsne$Y)
  row.names(points) =  row.names(map)
  colnames(points) = c("x", "y") #命名行名
  stress= NULL
}


#-----三种方法整体差异分析
# 三种方法计算群落差异分析函数
title1 = MicroTest(ps  = ps1_rela,Micromet = Micromet,dist = opts$distance)
title1

# 三种方法两两不比较

pairResult = pairMicroTest(ps = ps1_rela,Micromet = Micromet,dist = opts$distance)


#------------------------------------------plot
map = as.data.frame(sample_data(ps1_rela))
map$Group = as.factor(map$Group)
colbar <- length(levels(map$Group))

points = cbind(points, map[match(rownames(points), rownames(map)), ])
#write.table(points,"pcoa_bray_curtis.txt",quote = FALSE,row.names = F,
#           col.names = T,sep = "\t")
head(points)
points$ID = row.names(points)

mi = colorRampPalette(c( "#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
                         "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
                         "#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)
mi
# jpeg(file="./a2_bray_PCOA.jpeg")
if (method %in% c("DCA", "CCA", "RDA",  "MDS", "PCoA","PCA","LDA")) {
  p2 <-ggplot(points, aes(x=x, y=y, fill = Group)) +
    geom_point(alpha=.7, size=5, pch = 21) +
    labs(x=paste(method," 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste(method," 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),title=title1)+ #,title=title1
    stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))+
    #stat_ellipse( linetype = 1,level = 0.8)+
    #geom_text_repel(aes(label=points$id),size=4)+
    scale_colour_manual(values = mi,guide = guide_legend(title = NULL))+
    scale_fill_manual(values = mi,guide = guide_legend(title = NULL))+
    #labs(title = "toamto hea and dis")+
    guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
  p2
  # points$id=row.names(points)

  p2 = p2+theme_bw()+

    #scale_y_continuous(expand = c(0,0))+
    geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),

      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
      #legend.position = "none"#是否删除图例

    )
  p2
  library(ggrepel)
  p3 = p2+geom_text_repel( aes(label=points$ID),size=4)#?stat_ellipse

  p3

}




if (method %in% c("NMDS","t-sne")) {
  p2 <-ggplot(points, aes(x=x, y=y, fill = Group)) +
    geom_point(alpha=.7, size=5, pch = 21) +
    labs(x=paste(method,"1", sep=""),
         y=paste(method,"2",sep=""),
         title=stress)+
    stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))+
    #stat_ellipse( linetype = 1,level = 0.8)+
    #geom_text_repel(aes(label=points$id),size=4)+
    scale_colour_manual(values = mi,guide = guide_legend(title = NULL))+
    scale_fill_manual(values = mi,guide = guide_legend(title = NULL))+
    #labs(title = "toamto hea and dis")+
    guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
  p2
  # points$id=row.names(points)

  p2 = p2+theme_bw()+

    #scale_y_continuous(expand = c(0,0))+
    geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),

      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold")
      #legend.position = "none"#是否删除图例

    )
  library(ggrepel)
  p3 = p2+geom_text_repel( aes(label=points$ID),size=4)#?stat_ellipse

  p3
  if (method %in% c("t-sne")) {
   p2 = p2 +labs(x=paste(method,"1", sep=""),
            y=paste(method,"2",sep=""),
            title=title)
   p3 = p3 +labs(x=paste(method,"1", sep=""),
                 y=paste(method,"2",sep=""),
                 title=title)
  }
  p2

}
# ,pairResult,title1

return(list(p2,points,p3,pairResult,title1))

}





