# 生成STAMP输入物种分类级汇总表 format to stamp
#
# This is the format transform function named 'format2stamp'

#' @title Format otutab and taxonomy into each level for STAMP
#' @description Input otutab, taxonomy, and manual set abundance threshold.
#' dplyr merge taxonomy.
#' @param otutab raw or normalized OTU table
#' @param taxonomy Taxonomy include seven taxonomy level in tsv format
#' @param thre threshold of OTU relative abundance, default 1 read per Million
#' @param prefix prefix for output files
#' @details
#' By default, written 8 files
#' The files list as the following:
#' \itemize{
#' \item{1-7: Kingdom to Species together include 7 levels}
#' \item{8: filtered OTUs}
#' }
#' @return nothing, all table 1-8 saving in output directory.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#' Jingying Zhang, Yong-Xin Liu, et. al.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice
#' Nature Biotechnology. (2019) 37, 676-684.  DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#' @seealso format2lefse
#' @examples
#' # Set 4 parameters: 2 input files, and 2 parameters
#' format2stamp(otutab, taxonomy, thre = 0.01,  prefix = "tax")
#' @export

format2stamp <- function(otutab, taxonomy, thre = 0.0001, prefix = "tax") {

  # 内置数据测试函数
  # otutab=otutab
  # taxonomy=taxonomy
  # thre=0.0001
  # prefix="tax"

  # 数据交叉筛选
  idx = rownames(otutab) %in% rownames(taxonomy)
  otutab = otutab[idx,]
  tax = taxonomy[rownames(otutab),]

  # 依赖关系检测与安装
  p_list = c("dplyr")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 标准化，并筛选高丰度菌均值最小百万分之一0.0001%
  norm = t(t(otutab)/colSums(otutab,na=T))*100
  dim(norm)
  colSums(norm)
  idx = rowMeans(norm) > thre
  HA = norm[idx,]
  dim(HA)
  colSums(HA)
  # 保在筛选后的OTU表
  write.table(paste("OTU\t",  sep=""), file=paste(prefix, "_8OTU", thre, ".txt", sep = ""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
  suppressWarnings(write.table(HA, file=paste(prefix, "_8OTU", thre, ".txt", sep = ""), append = T, sep="\t", quote=F, row.names=T, col.names=T))
  # 数据筛选并排序，要求每个OTU必须有注释，可以为空
  tax = tax[rownames(HA),]

  # 转换为等级|连接格式
  tax$Phylum=paste(tax$Kingdom,tax$Phylum,sep = "|")
  tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
  tax$Order=paste(tax$Class,tax$Order,sep = "|")
  tax$Family=paste(tax$Order,tax$Family,sep = "|")
  tax$Genus=paste(tax$Family,tax$Genus,sep = "|")
  tax$Species=paste(tax$Genus,tax$Species,sep = "|")
  # head(tax)

  # 按Kingdom合并
  grp <- tax[rownames(tax), "Kingdom", drop=F]
  merge=cbind(HA, grp)
  HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
  colnames(HA_Kingdom)[1]="Kingdom"
  # 添加3有效小更简洁和节约空间 round(matrix, digits = 3)
  write.table(HA_Kingdom, file=paste(prefix, "_1Kingdom.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Kingdom)[1]="Class"

  # 按Phylum合并
  grp <- tax[rownames(tax), "Phylum", drop=F]
  merge=cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  colnames(HA_Phylum)[1]="Phylum"
  write.table(HA_Phylum, file=paste(prefix, "_2Phylum.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Phylum)[1]="Class"

  # 按Class合并
  grp <- tax[rownames(tax), "Class", drop=F]
  merge=cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  colnames(HA_Class)[1]="Class"
  write.table(HA_Class, file=paste(prefix, "_3Class.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Class)[1]="Class"

  # 按Order合并
  grp <- tax[rownames(tax), "Order", drop=F]
  merge=cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  colnames(HA_Order)[1]="Order"
  write.table(HA_Order, file=paste(prefix, "_4Order.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Order)[1]="Class"

  # 按Family合并
  grp <- tax[rownames(tax), "Family", drop=F]
  merge=cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  colnames(HA_Family)[1]="Family"
  write.table(HA_Family, file=paste(prefix, "_5Family.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Family)[1]="Class"

  # 按Genus合并
  grp <- tax[rownames(tax), "Genus", drop=F]
  merge=cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  colnames(HA_Genus)[1]="Genus"
  write.table(HA_Genus, file=paste(prefix, "_6Genus.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Genus)[1]="Class"

  # 按Species合并
  grp <- tax[rownames(tax), "Species", drop=F]
  merge=cbind(HA, grp)
  HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
  colnames(HA_Species)[1]="Species"
  write.table(HA_Species, file=paste(prefix, "_7Species.txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
  colnames(HA_Species)[1]="Class"
}
