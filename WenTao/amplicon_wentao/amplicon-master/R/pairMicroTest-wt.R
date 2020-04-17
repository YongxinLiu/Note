
# 群落两两比差异比较
# result = pairMicroTest(ps = ps,dist = "bray",Micromet = "MRPP")

# ps = readRDS("./ps_liu.rds")
# unif <- phyloseq::distance(ps, method="bray")
# unif
# ps = ps
# dist = "bray"
# Micromet = "MRPP"
# unif <- phyloseq::distance(ps, method="bray")
# 经过调整后，适合了shell，但是该函数在R中的使用就不是很快了，需要将otu等表格转化phyloseq，其次要计算距离，然后才能使用
# unif <- phyloseq::distance(ps, method="bray")

# pairMicroTest (ps = ps,Micromet = "anosim",dist = 8)

pairMicroTest = function(ps = ps,Micromet = "anosim",dist = 8){
  library(phyloseq)
  dist_methods <- unlist(phyloseq::distanceMethodList)
  dist =   dist_methods[dist]
  ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  library(vegan)
  #-准备矩阵和分组文件
  map = as.data.frame(sample_data(ps1_rela))
  
  aa = levels(map$Group)
  aa
  aaa = combn(aa,2)
  aaa 
  dim(aaa)[2]
  
  # 构建三个空列
  ID = rep("a",dim(aaa)[2])
  R = rep("a",dim(aaa)[2])
  P = rep("a",dim(aaa)[2])
  # i = 2
  # as = as.matrix(unif)
  for (i in 1:dim(aaa)[2]) {
    print(i)
    Desep_group = aaa[,i]
    map = as.data.frame(sample_data(ps1_rela))
    head(map)
    map$ID = row.names(map)
    maps<- dplyr::filter(map,Group %in% Desep_group)
    row.names(maps) = maps$ID
    ps_sub = ps1_rela
    sample_data( ps_sub ) = maps
    # ps_sub <- phyloseq::subset_samples(ps1_rela,Group %in% Desep_group);ps_sub 
    
    
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE);ps_sub
    map = as.data.frame(sample_data(ps_sub))
    unif <- phyloseq::distance(ps_sub, method=dist)
    # # 取子集
    # 
    # unif = as[map$ID,map$ID]
    # 
    
    # print(ps_sub)
    
    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map$Group) 
      as1 = round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as1, sep = "")
      # R[i] = R2
      R2
      p_v = paste("p: ",round(mrpp$Pvalue,3), sep = "")
      p_v
      # P[i] = p_v
      
    }
    
    if (Micromet == "anosim") {
      dat.ano = anosim(unif, map$Group) 
      a = round(dat.ano$statistic,3)
      R2 <- paste("ANOSIM.r ",a, sep = "")
      R[i] = R2
      p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
      P[i] = p_v
      # title = paste(R2," ",p_v, sep = "")
      # title
      
    }
  gg =map$Group

    if (Micromet == "adonis") {
     
      ado =  adonis(unif~gg,permutations = 999)

      a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
      R2 <- paste("adonis:R ",a, sep = "")
      R[i] = R2
      b = as.data.frame(ado$aov.tab[6])[1,1]
      p_v = paste("p: ",b, sep = "")
      P[i] = p_v
      # title = paste(R2," ",p_v, sep = "")
      # title
      # print(i)
      # print(R)
    }
  print("iii")
    ID[i] = paste(Desep_group[1],Desep_group[2],sep = "_VS_")
    P[i] = p_v
    R[i] = R2
  }
  P
  R
  result = data.frame(ID = ID,stat = R,p = P)
  head(result)
  
  return(result)
}


