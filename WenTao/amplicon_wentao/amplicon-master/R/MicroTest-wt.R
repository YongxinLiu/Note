

# 经过调整后，适合了shell，但是该函数在R中的使用就不是很快了，需要将otu等表格转化phyloseq，其次要计算距离，然后才能使用
# unif <- phyloseq::distance(ps, method="bray")

# ps  = ps1_rela
# Micromet = Micromet
# dist = dist

MicroTest = function(otu = NULL,map = NULL,ps = NULL,group = "Group",Micromet = "MRPP",dist = 8){
  
  library(phyloseq)
  dist_methods <- unlist(phyloseq::distanceMethodList)
  dist =   dist_methods[dist]

  # # 数据导入
  # if (is.null(otu)&is.null(map)) {
  #   ps = ps
  # }else {
  #   #导入otu表格
  #   otu = as.matrix(otu)
  #   str(otu)
  #   # #导入注释文件
  #   # tax = as.matrix(tax)
  #   # taxa_names(tax)
  # 
  #   #导入分组文件
  #   colnames(map) = gsub(group,"AA", colnames(map))
  # 
  #   map$Group = map$AA
  #   map$Group = as.factor(map$Group )
  #   map$Group
  #   # #导入进化树
  #   # tree = read.tree("./otus.tree")
  #   # tree
  #   # ?phyloseq
  #   # ps <- phyloseq(
  #   #   otu_table(otu, taxa_are_rows=TRUE),
  #   #                sample_data(map)
  #   #                # phy_tree(tree)
  #   # )
  #   #
  # 
  #   otu = phyloseq::otu_table(otu,taxa_are_rows=TRUE)
  #   map = phyloseq::sample_data(map)
  #   ps = phyloseq::merge_phyloseq(otu,map)
  # ps
  # }
  # 
  # # dist_methods <- unlist(distanceMethodList)
  # if (dist %in% c("unifrac" , "wunifrac",  "dpcoa")) {
  #   phy_tree(ps) = tree
  # }
  # 
  # ps = readRDS("./ps_liu.rds")
  
  # 求取相对丰度
  ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela 

  library(vegan)
  # -准备矩阵和分组文件
  
  map = as.data.frame(sample_data(ps1_rela))
  # ?distance
  # unif<- distance(ps1_rela , method=dist)
  # unif = vegdist(t(otu), methodC="bray")
  # unif <- distance(otu, method=dist)
  
  
  unif <- phyloseq::distance(ps, method=dist)

  if (Micromet == "MRPP") {
    dat.mrpp = mrpp(unif, map$Group)
    a = round(dat.mrpp$delta,3)
    R2 <- paste("MRPP.delta ",a, sep = "")
    p_v = paste("p: ",round(dat.mrpp$Pvalue,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

  if (Micromet == "anosim") {
    dat.ano = anosim(unif, map$Group)
    a = round(dat.ano$statistic,3)
    R2 <- paste("ANOSIM.r ",a, sep = "")
    p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1

  }
  if (Micromet == "adonis") {
    ado =  vegan::adonis(unif  ~ map$Group,permutations = 999)
    a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 <- paste("adonis:R ",a, sep = "")
    b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",b, sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

  return(title1)
}
