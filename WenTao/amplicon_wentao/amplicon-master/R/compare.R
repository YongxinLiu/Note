# 差异比较 Difference comparison
#
# This is the function named 'compare'
# which output difference comparison table by wilcoxon, t.test or edgeR
#
#' @title Difference comparison
#' @description Input feature table, group info and comparing pair
#' Features with relative abundance, p-value and FDR for pair.
#' @param data data matrix
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param group column name for groupID.
#' @param compare_pair paired groups linked with dash.
#' @param RA threshold for relative abundance.
#' @param method method of test.
#' @param pvalue threshold of pvalue.
#' @param fdr threshold of fdr.
#' @details
#' By default, return table
#' The available test methods include the following:
#' \itemize{
#' \item{most used methods: wilcoxon test, t.test}
#' \item{other used indices: edgeR, DESeq2}
#' }
#' @return Difference comparison table.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, et. al. (2019).
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology 37, 676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso compare_vocalno compare_heatmap compare_manhattan
#' @examples
#' # Set 4 parameters
#' compare(data, metadata, group, compare, method)
#' # Set 4 parameters
#' compare(data, metadata, group, compare, method, RA, pvalue, fdr)
#' @export



compare <- function(data = otutab, metadata = metadata, group = "genotype", compare_pair = "KO-WT", method = "wilcox", RA = 0.01, pvalue = 0.05, fdr = 0.1) {

# 依赖关系检测与安装
# p_list = c("edgeR", "DESeq2")
# for(p in p_list){
#   if (!requireNamespace(p)){
#   install.packages(p)}
#   library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

#----测试默认参数#----
# data = otutab
# metadata = metadata
# group = "genotype"
# method = "edgeR"
# compare_pair = "KO-WT"
# RA = 0.05
# pvalue = 0.05
# fdr = 0.1
#----测试命令行参数#----
# data = otutab
# metadata = metadata
# group = opts$group
# compare_pair = opts$compare
# method = opts$method
# RA = opts$threshold
# pvalue = opts$pvalue
# fdr = opts$fdr
#----交叉筛选#----
idx = rownames(metadata) %in% colnames(data)
metadata = metadata[idx,]
data = data[, rownames(metadata)]

#----丰度过滤#----
# 标准化为百分比
norm = t(t(data)/colSums(data,na=T)*100)
# 按丰度筛选标准化特征表和原始值
idx = rowMeans(norm) > RA
norm = norm[idx, ]
colSums(norm)
data = data[idx, ]

#----差异比较组筛选#----
group_list = strsplit(compare_pair,'-')[[1]]
metadata$group=metadata[,group]
idx = metadata$group %in% group_list
sub_metadata = metadata[idx,]
sub_dat=as.matrix(data[, rownames(sub_metadata)])


#----edgeR差异比较#----
if (method == "edgeR"){
  print("Your are using edgeR test!")
  library(edgeR)
  d = DGEList(counts=sub_dat,group=factor(sub_metadata$group))
  d = calcNormFactors(d)
  # check samples is in right groups
  # d$samples
  # metadata.mat = model.matrix(~factor(sub_metadata$group))
  metadata.mat = model.matrix(~ 0 + factor(sub_metadata$group))
  rownames(metadata.mat) = colnames(sub_dat)
  colnames(metadata.mat) = levels(factor(sub_metadata$group))
  DAF = estimateDisp(d,metadata.mat)
  fit = glmFit(DAF,metadata.mat)
  BvsA <- makeContrasts(contrasts = compare_pair, levels = metadata.mat)
  lrt = glmLRT(fit,contrast=BvsA)
  nrDAF=as.data.frame(topTags(lrt, n=nrow(sub_dat)))
  # nrDAF=as.data.frame(nrDAF)
  nrDAF = nrDAF[, -3]
}else if (method == "wilcox"){
  print("Your are using Wilcoxon test!")
  idx = sub_metadata$group %in% group_list[1]
  GroupA = norm[,rownames(sub_metadata[idx,,drop=F])]
  idx = sub_metadata$group %in% group_list[2]
  GroupB = norm[,rownames(sub_metadata[idx,,drop=F])]

  nrDAF = data.frame(list=rownames(norm), row.names =rownames(norm) )
  # 对每行Feature进行秩合检验
  for ( i in 1:dim(nrDAF)[1]){
    FC = (mean(GroupA[i,])+0.0000001)/(mean(GroupB[i,])+0.0000001)
    nrDAF[i,2]=log2(FC)
    nrDAF[i,3]=log2(max(c(GroupA[i,],GroupB[i,]))*10000)
    nrDAF[i,4]= suppressWarnings(wilcox.test(as.numeric(GroupA[i,]),as.numeric(GroupB[i,]))$p.value)
  }
  nrDAF=nrDAF[,-1]
  colnames(nrDAF)=c("logFC", "logCPM", "PValue")
  nrDAF$FDR = p.adjust(nrDAF$PValue, method="fdr", dim(nrDAF)[1])
}else if (method == "t.test"){
  print("Your are using T test!")
  idx = sub_metadata$group %in% group_list[1]
  GroupA = norm[,rownames(sub_metadata[idx,,drop=F])]
  idx = sub_metadata$group %in% group_list[2]
  GroupB = norm[,rownames(sub_metadata[idx,,drop=F])]

  nrDAF = data.frame(list=rownames(norm), row.names =rownames(norm) )
  # 对每行Feature进行秩合检验
  for ( i in 1:dim(nrDAF)[1]){
    FC = (mean(GroupA[i,])+0.0000001)/(mean(GroupB[i,])+0.0000001)
    nrDAF[i,2]=log2(FC)
    nrDAF[i,3]=log2(max(c(GroupA[i,],GroupB[i,]))*10000)
    nrDAF[i,4]= suppressWarnings(t.test(as.numeric(GroupA[i,]),as.numeric(GroupB[i,]))$p.value)
  }
  nrDAF=nrDAF[,-1]
  colnames(nrDAF)=c("logFC", "logCPM", "PValue")
  nrDAF$FDR = p.adjust(nrDAF$PValue, method="fdr", dim(nrDAF)[1])
  }else{
  print("Unknown method!")
}


# 结果调整3位数，并分三类
nrDAF$logFC=round(nrDAF$logFC,3)
nrDAF$logCPM=round(nrDAF$logCPM,3)
# nrDAF$PValue=round(nrDAF$PValue,5)
# nrDAF$FDR=round(nrDAF$FDR,5)
nrDAF$level = ifelse(nrDAF$logFC>0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Enriched",
                     ifelse(nrDAF$logFC<0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Depleted",
                            "NotSig"))
# nrDAF$level=factor(nrDAF$level,levels = c("Enriched","Depleted","NotSig"))
print(table(nrDAF$level))
# head(nrDAF)

#----添加均值#----
# Add MeanA and MeanB in percentage
# calculate groupA mean
A_list = subset(sub_metadata, group %in% group_list[1])
A_norm = norm[, rownames(A_list)]
A_mean = as.data.frame(rowMeans(A_norm))
colnames(A_mean)=c("MeanA")
# calculate groupB mean
B_list = subset(sub_metadata, group %in% group_list[2])
B_norm = norm[, rownames(B_list)]
B_mean = as.data.frame(rowMeans(B_norm))
colnames(B_mean)=c("MeanB")
# merge and reorder
Mean = round(cbind(A_mean, B_mean, A_norm, B_norm), 3) #
Mean = Mean[rownames(nrDAF),]

# 修正列名
colnames(nrDAF)[1] = "log2FC"
colnames(nrDAF)[2] = "log2CPM"
output=cbind(nrDAF,Mean)

}
