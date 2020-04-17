# 生成LEfSe输入分类级汇总表 format to LEfSe
#
# This is the first function named 'format2lefse'

#' @title Format otutab, metadata and taxonomy into LEfSe format with level prefix p_
#' @description Input otutab, taxonomy and metadata, and manual set abundance threshold and metadata column names.
#' dplyr merge taxonomy.
#' @param otutab OTUID in row and normalized OTU table
#' @param taxonomy Taxonomy include seven taxonomy level in tsv format
#' @param metadata matrix or dataframe, including sampleID and groupID
#' @param thre threshold of OTU abundance, especially for lefse create properly cladogram
#' @param groupID column name for group ID.
#' @details
#' By default, written 2 file
#' @return LEfSe input with or without taxonomy letter.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#' Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).  DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#' @seealso format2stamp
#' @examples
#' # Set six parameters: three input files, and three parameters
#' format2lefse(otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", output = "LEfSe.txt")
#' @export

format2lefse2 <- function(otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", output = "LEfSe.txt") {

  # 内置数据测试函数
  # otutab=otutab
  # taxonomy=taxonomy
  # metadata=metadata
  # thre=0.01
  # groupID = "genotype"
  # output = "LEfSe.txt"

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
  # dim(norm)
  # colSums(norm)
  idx = rowMeans(norm) > thre
  HA = norm[idx,]
  # dim(HA)
  # colSums(HA)
  # 数据筛选并排序，要求每个OTU必须的注释，可以为空
  tax = tax[rownames(HA),]

  # 级别添加标签
  tax$Kingdom = paste0("k__",tax$Kingdom)
  tax$Phylum = paste0("p__",tax$Phylum)
  tax$Class = paste0("c__",tax$Class)
  tax$Order = paste0("o__",tax$Order)
  tax$Family = paste0("f__",tax$Family)
  tax$Genus = paste0("g__",tax$Genus)

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
  colnames(HA_Kingdom)[1]="Class"

  # 按Phylum合并
  grp <- tax[rownames(tax), "Phylum", drop=F]
  merge=cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  colnames(HA_Phylum)[1]="Class"

  # 按Class合并
  grp <- tax[rownames(tax), "Class", drop=F]
  merge=cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  colnames(HA_Class)[1]="Class"

  # 按Order合并
  grp <- tax[rownames(tax), "Order", drop=F]
  merge=cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  colnames(HA_Order)[1]="Class"

  # 按Family合并
  grp <- tax[rownames(tax), "Family", drop=F]
  merge=cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  colnames(HA_Family)[1]="Class"

  # 按Genus合并
  grp <- tax[rownames(tax), "Genus", drop=F]
  merge=cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  colnames(HA_Genus)[1]="Class"

  # 合并6个分类级
  all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)

  # 将选定的分组列统一命名为group
  metadata$group = metadata[, groupID]

  # 修改样品名为组名，用实验设计中的group列替换样品名
  colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]],]$group)

  # 保存结果为lefse
  write.table(all, file=paste(output, sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
}
