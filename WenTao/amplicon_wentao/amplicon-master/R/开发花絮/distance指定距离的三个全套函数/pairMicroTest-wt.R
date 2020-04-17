
# 群落两两比差异比较
# result = pairMicroTest(ps = ps,dist = "bray",Micromet = "MRPP")



pairMicroTest = function(ps = ps,dist = "bray",Micromet = "anosim"){
  
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
    unif <- phyloseq::distance(ps_sub  , method=dist, type="samples")
    # print(ps_sub)
    
    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map$Group) 
      as = round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as, sep = "")
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
    if (Micromet == "adonis") {
      ado =  vegan::adonis(unif  ~ map$Group,permutations = 999)
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


