
# #清空内存
# rm(list=ls()) 
# library("tidyverse")
# # ps
# # allp = c("Phylum","Class","Order","Family","Genus")
# # ## 首先按照属进行合并
# # taxGlomRank = "Genus"
# # taxGlomRank = "Phylum"
# # N = 30
# ## 导入数据
# ps = readRDS("./a3_DADA2_table//ps.rds")
# ps
# result = BubblePlot(ps = ps ,N = 30,taxGlomRank = "Phylum",group = "SampleType")
# 
# # 提取气泡图:样本为纵向
# p1 = result[[1]]
# # 样本为横向
# p2 = result[[2]]
# p2
# 
# # 提取作图数据
# 
# dataplot = result[[3]]
# dataplot
# ps = ps 
# N = 30
# taxGlomRank = "Phylum"
# group = "Group"
# 导入otu表格
# otu = read.delim("./data/otutab.txt",row.names = 1)
# #导入注释文件
# tax = read.delim("./data/taxonomy.txt",row.names = 1)
# head(tax)
# #导入分组文件
# map = read.delim("./data/metadata.tsv",row.names = 1)
# head(map)
# 
# 
# result = BubblePlot(otu = otu,tax = tax,map = map ,N = 30,taxGlomRank = "Phylum",group = "Group")
# 
# # 提取气泡图:样本为纵向
# p1 = result[[1]]
# # 样本为横向
# p2 = result[[2]]
# p2


# otu = otu
# tax = tax
# map = map
# ps = NA
# N = 30
# taxGlomRank = "Phylum"
# group = "Group"
# axis_ord = NA

BubblePlot = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,axis_ord = NULL,N = 30,taxGlomRank = "Phylum",group = "Group"){
  library(phyloseq)
  library(tidyverse)
  # 功能函数
  ### 提取OTU表格
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  ### 添加OTU注释信息
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)
    
    return(as(tax,"matrix"))
  }
  
  
  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }
  
  if (is.null(otu)&is.null(tax)&is.null(map)) {
    ps = ps
    
  }
  
  if (!is.null(otu)|!is.null(tax)|!is.null(map) ) {
   
    head(otu)
    otu = as.matrix(otu)
    str(otu)
    
    tax = as.matrix(tax)
    # taxa_names(tax)
    
    
    colnames(map) = gsub(group,"AA", colnames(map))
    
    map$Group = map$AA
    map$Group = as.factor(map$Group )
    map$Group 
    # #导入进化树
    # tree = read.tree("./otus.tree")
    # tree
    
    ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE), 
                   sample_data(map) ,
                   tax_table(tax)
                   # phy_tree(tree)
    )
  }
  ps
  map = as.data.frame(sample_data(ps))
  head(map)
  colnames(map) = gsub(group,"AA", colnames(map))
  
  map$Group = map$AA
  map$Group = as.factor(map$Group )
  sample_data(ps) = map
  
  
  ps_sub = tax_glom(ps, taxrank = taxGlomRank)
  
  #----------------按照分组，求取每组的平均丰度
  otu_table = as.data.frame(t(vegan_otu(ps_sub)))
  head(otu_table)
  design = as.data.frame(sample_data(ps_sub))
  ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
  OTU = as.matrix(otu_table)
  norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
  norma = norm %>% 
    t() %>% as.data.frame()
  #数据分组计算平均值
  iris.split <- split(norma,as.factor(design$Group))
  
  iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
  # 组合结果
  norm2 <- do.call(rbind,iris.apply)%>% 
    t() 
  norm2 = as.data.frame(norm2)
  norm2$mean=apply(norm2,1,mean)
  norm2$ID = row.names(norm2)
  colnames(norm2)
  ##按照mean、列进行排序desc设置从大到小排序
  norm3<- arrange(norm2, desc(mean))
  row.names(norm3) = norm3$ID
  norm3$ID = NULL
  
  
  
  #---------------提取展示的物种
  ### 提取前30个属
  wt = head(norm3,n = N)
  wt$mean = NULL
  
  
  tax_table = as.data.frame(vegan_tax(ps_sub))
  head(tax_table)
  wt
  wt_tax = merge(wt,tax_table,by = "row.names",all = F)
  head(wt_tax)
  row.names(wt_tax) = wt_tax[,taxGlomRank]
  wt = wt_tax[,c(colnames(wt))]
  head(wt)
  wt$ID = row.names(wt) 
  library(reshape2)
  pcm = melt(wt, id = c("ID"))
  head(pcm)
  
  pcm$ID <- factor(pcm$ID,levels=unique(pcm$ID))
  
  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  #----样本在y轴上
  p1 = ggplot(pcm, aes(x = ID, y = variable)) + 
    geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,17), breaks = c(1,10,50,75)) + 
    labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
    theme(legend.key=element_blank(), 
          axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
          axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
          legend.text = element_text(size = 10, face ="bold", colour ="black"), 
          legend.title = element_text(size = 12, face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.position = "right") +  
    scale_fill_manual(values = colours, guide = FALSE) + 
    scale_y_discrete(limits = rev(levels(pcm$variable))) 
  
  p1
  
  #---样本在x轴上
  p2 = ggplot(pcm, aes(y = ID, x = variable)) + 
    geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,17), breaks = c(1,10,50,75)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    theme(legend.key=element_blank(), 
          axis.text.y = element_text(colour = "black", size = 12, face = "bold", angle = 0, vjust = 0.3, hjust = 1), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 11), 
          legend.text = element_text(size = 10, face ="bold", colour ="black"), 
          legend.title = element_text(size = 12, face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
          legend.position = "right") +  
    scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable))) 
  
  p2
  
  return(list(p1,p2,pcm))
}








