

#---------纯R语言绘制稀释曲线

#-调用ps对象
library(tidyverse)
ps = readRDS("../ps_liu.rds")


#------------构造随机抽样函数
vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
otb = as.data.frame(t(vegan_otu(ps)))



#--------------开始抽平操作；
#设置抽平数量
N = 100
#-
i = 1
ii =1
a = c()
aa = c()
for (i in 1:dim(otb)[1]) {
  
  b= paste(row.names(otb[i,]),1:otb[i,1],sep = "&")
  a = c(a,b)
}

ss = sample(a,N,replace = FALSE)
str_split(ss, pattern = "&", n = Inf, simplify = FALSE)

sample.names <- sapply(str_split(ss, pattern = "&", n = Inf, simplify = FALSE), `[`, 1)


aa = as.data.frame(table(sample.names))
sum(aa$Freq)
colnames(aa)= c("ID",colnames(otb)[ii])
head(aa)
row.names(aa) = aa$ID
aa$ID = NULL


for (ii in 2:dim(otb)[2]) {
  
  for (i in 1:dim(otb)[1]) {
    b= paste(row.names(otb[i,]),1:otb[i,ii],sep = "&")
    a = c(a,b)
  }
  
  ss = sample(a,N,replace = FALSE)
  # str_split(ss, pattern = "&", n = Inf, simplify = FALSE)
  
  sample.names <- sapply(str_split(ss, pattern = "&", n = Inf, simplify = FALSE), `[`, 1)
  bb = as.data.frame(table(sample.names))
  # bb[2]
  colnames(bb)= c("ID",colnames(otb)[ii])
  head(bb)
  row.names(bb) = bb$ID
  bb$ID = NULL
  head(bb)
  aa = merge(aa,bb,by = "row.names",all = TRUE)
  row.names(aa) = aa$Row.names
  aa$Row.names = NULL
  
}


dim(aa)
aa[is.na(aa)] = 0

head(aa)



#--------------结束抽平--------------------------------------------------------------------------------

