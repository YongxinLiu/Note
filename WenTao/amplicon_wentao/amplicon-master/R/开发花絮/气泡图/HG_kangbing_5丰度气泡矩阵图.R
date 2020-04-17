#########相对丰度气泡图




#-----------导入数据
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)

L2 = read.table("otu_table_tax_L2.txt", header=T,  sep="\t") 
head(L2)


# 过滤数据并排序，只有定义为行名是才可以排序
rownames(design)=design$SampleID2
idx = rownames(design) %in% colnames(L2) 
idx
sub_design = design[idx,]
count = L2[, rownames(sub_design)]
head(count)









library(limma)
#下面来筛选差异otu
design.mat = model.matrix(~ 0 + sub_design$SampleType)
colnames(design.mat)=levels(design$SampleType)
#可以同时设置好几组比较
contrast.matrix <- makeContrasts(CSF-CRF, levels=design.mat)
#行线性模型拟合
fit <- lmFit(count, design.mat)
#根据对比模型进行差值计算T-test对数据进行计算
fit2 <- contrasts.fit(fit, contrast.matrix)
#贝叶斯检验
fit2 <- eBayes(fit2)

results<-decideTests(fit2, method="global", adjust.method="BH", p.value=0.05, lfc=0)
summary(results)
x<-topTable(fit2, coef=1, number=10000, adjust.method="BH", sort.by="p", resort.by=NULL)
head(x)
x$levelLF = as.factor(ifelse(x$adj.P.Val < 0.05 & x$logFC > 0, "enriched",ifelse(x$adj.P.Val < 0.05 & x$logFC < 0, "nosig","nosig")))
x$levelB80 = as.factor(ifelse(x$adj.P.Val < 0.05 & x$logFC > 0, "nosig",ifelse(x$adj.P.Val < 0.05 & x$logFC < 0, "depleted","nosig")))


#######计算相对丰度均值
# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)
normB80=norm[1:6]
head(normB80)
normB80$meanB80=apply(normB80,1,mean)
###
normLF=norm[7:12]
head(normLF)
normB80$meanLF=apply(normLF,1,mean)
head(normB80)
normB80[grep(".fq|Row.names",colnames(normB80))]<-NULL

index = merge(normB80,x, by="row.names",all=F)
head(index)

index2=data.frame(name=index$Row.names,LF=index$meanLF,B80=index$meanB80)
head(index2)
index3=data.frame(name=index$Row.names,LF=index$levelLF,B80=index$levelB80)
head(index3)

###########
library (ggplot2)
library (reshape2)
## 利用reshape2将数据框从宽型重构为长型

tax <- melt (index3, id="name")
head(tax)
colnames(tax)=c("name","break1","fengzu")

#########
fengdu <- melt (index2, id="name")
head(fengdu)
colnames(fengdu)=c("name","break1","fengdu")
#########


#########
## 利用ggplot2的散点图作图
## 样品品映射为x轴，属名映射为y轴
## 丰度映射为气泡大小
######将数据转化#wt2<-sqrt(wt)
fengdu$log10=-log10(fengdu$fengdu+0.000001)
head(fengdu)
fengdu$fengzu=tax$fengzu
#注意必须转化为因子
fengdu$break1=factor(fengdu$break1)
fengdu$name=factor(fengdu$name)
#####position = position_dodge(0)设置倾斜度yintercept = 10,xintercept = 10#, group=tax$fengzu
mi=c("red","green","#FFFFB3")

p <- ggplot (fengdu, aes (x=break1, y=name,size=log10,fill=fengzu))+ 
  geom_point(shape=21,colour="black" )+scale_size_area(max_size = 3)+scale_fill_manual(values =mi)+
  #scale_fill_gradient2(low = "red", high = "blue")+
  #geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1.5,colour="white")+
  geom_hline(data=fengdu,aes(yintercept=1.5:94.5),colour="white")

p
p +theme(axis.text.x =element_text(angle = 90,vjust=-0.05),
         axis.text.y =element_text(size=6), 
         panel.background = element_rect(fill = "grey90"),
)+
  coord_fixed(ratio = 1.2)




########
##########纲水平相对丰度气泡图##################
##########纲水平相对丰度气泡图##################
##########纲水平相对丰度气泡图##################
##########纲水平相对丰度气泡图##################


setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli")
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)

setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/taxa_summary")
L3 = read.table("otu_table_tax_L3.txt", header=T, row.names= 1, sep="\t") 
head(L3)


# 过滤数据并排序，只有定义为行名是才可以排序
rownames(design)=design$SampleID2
idx = rownames(design) %in% colnames(L3) 
idx
sub_design = design[idx,]
count = L3[, rownames(sub_design)]
head(count)

#下面来筛选差异otu
design.mat = model.matrix(~ 0 + sub_design$SampleType)
colnames(design.mat)=levels(design$SampleType)
#可以同时设置好几组比较
contrast.matrix <- makeContrasts(LF-B80, levels=design.mat)
#行线性模型拟合
fit <- lmFit(count, design.mat)
#根据对比模型进行差值计算T-test对数据进行计算
fit2 <- contrasts.fit(fit, contrast.matrix)
#贝叶斯检验
fit2 <- eBayes(fit2)
#"fdr" (equivalent to "BH");lfc(log2-fold-change )
results<-decideTests(fit2, method="global", adjust.method="BH", p.value=0.05, lfc=0)
summary(results)
x<-topTable(fit2, coef=1, number=10000, adjust.method="BH", sort.by="p", resort.by=NULL)
head(x)
x$levelLF = as.factor(ifelse(x$adj.P.Val < 0.05 & x$logFC > 0, "enriched",ifelse(x$adj.P.Val < 0.05 & x$logFC < 0, "nosig","nosig")))
x$levelB80 = as.factor(ifelse(x$adj.P.Val < 0.05 & x$logFC > 0, "nosig",ifelse(x$adj.P.Val < 0.05 & x$logFC < 0, "depleted","nosig")))

#######计算相对丰度均值
# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)
normB80=norm[1:6]
head(normB80)
normB80$meanB80=apply(normB80,1,mean)
###
normLF=norm[7:12]
normB80$meanLF=apply(normLF,1,mean)
normB80[grep(".fq|Row.names",colnames(normB80))]<-NULL
head(normB80)
index = merge(normB80,x, by="row.names",all=F)
head(index)

index2=data.frame(name=index$Row.names,LF=index$meanLF,B80=index$meanB80)
head(index2)
index3=data.frame(name=index$Row.names,LF=index$levelLF,B80=index$levelB80)
head(index3)

###########
library (ggplot2)
library (reshape2)
## 利用reshape2将数据框从宽型重构为长型

tax <- melt (index3, id="name")
head(tax)
colnames(tax)=c("name","break1","fengzu")

#########
fengdu <- melt (index2, id="name")
head(fengdu)
colnames(fengdu)=c("name","break1","fengdu")
#########


#########
## 利用ggplot2的散点图作图
## 样品品映射为x轴，属名映射为y轴
## 丰度映射为气泡大小
######将数据转化#wt2<-sqrt(wt)
fengdu$log10=-log10(fengdu$fengdu+0.000001)
head(fengdu)
fengdu$fengzu=tax$fengzu
#注意必须转化为因子
fengdu$break1=factor(fengdu$break1)
fengdu$name=factor(fengdu$name)
#####position = position_dodge(0)设置倾斜度yintercept = 10,xintercept = 10#, group=tax$fengzu
mi=c("red","green","#FFFFB3")
library(ggthemes)
pdf("L2.pdf",width =2.8,height = 8.5) 
p <- ggplot (fengdu, aes (x=break1, y=name,size=log10,fill=fengzu))+ 
  geom_point(shape=21,colour="black" )+scale_size_area(max_size = 3)+scale_fill_manual(values =mi)+
  #scale_fill_gradient2(low = "red", high = "blue")+
  #geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1.5,colour="white")+
  geom_hline(data=fengdu,aes(yintercept=1.5:242.5),colour="white")
p

p + theme_few()+theme(axis.text.x =element_text(angle = 90,vjust=-0.05),
                      axis.text.y =element_text(size=6), 
                      panel.background = element_rect(fill = "grey90"),
)+
  coord_fixed(ratio = 0.6)

dev.off() 


