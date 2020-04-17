# 测试
# title1 = MicroTest(ps  = ps1_rela,Micromet = Micromet)
# title1

MicroTest = function(ps,Micromet = "MRPP",dist = dist){
  ps = ps
  ps1_rela  = transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  library(vegan)
  #-准备矩阵和分组文件
  map = as.data.frame(sample_data(ps1_rela))
  unif <- phyloseq::distance(ps1_rela , method=dist, type="samples")
  
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
