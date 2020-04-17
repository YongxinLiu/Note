# 绘制beta多样性组间统计
#
# This is the function named 'beta_pcoa_stat'
# which draw PCoA scatter plot with stat ellipse, and reture a ggplot2 object
#
#' @title Plotting beta diversity scatter plot
#' @description Input alpha index and metadata, and manual set alpha index and metadata column names.
#' ggplot2 show PCoA with color and stat ellipse.
#' @param dis_mat distance matrix, typical output of usearch -beta_div,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @param pairwise group pairwise compare, only groups less than 5, T or F.
#' @param pairwise_list a file include pairwise list.
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#' @seealso beta_cpcoa
#' @examples
#' # Set 3 parameters: dis_mat, metadata and groupID
#' beta_pcoa_stat(dis_mat = beta_bray_curtis, metadata, "genotype")
#' # Set 5 parameters: dis_mat, metadata, and groupID, pairwise FLASE, can set pairwise_list in file
#' beta_pcoa_stat(dis_mat, metadata, groupID = "genotype", pairwise = TRUE, pairwise_list = "doc/compare.txt")
#' @export

#两两adonis，
# beta_pcoa_stat(dis_mat = beta_bray_curtis, metadata, "genotype")


beta_pcoa_stat <- function(dis_mat, metadata, groupID = "genotype", pairwise = TRUE, pairwise_list = "doc/compare.txt") {

  # 依赖关系检测与安装
  p_list = c("vegan")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  metadata$group = metadata[[groupID]]

  # Compare each group beta by vegan adonis in bray_curtis
  da_adonis = function(sampleV){
    sampleA = as.matrix(sampleV$sampA)
    sampleB = as.matrix(sampleV$sampB)
    design2 = subset(metadata, group %in% c(sampleA,sampleB))
    if (length(unique(design2$group))>1) {
      sub_dis_table = dis_table[rownames(design2),rownames(design2)]
      sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
      adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000)
      adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
      # print(paste("In ",opts$type," pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
      adonis_pvalue = paste(sampleA, sampleB, adonis_pvalue, sep="\t")
      write.table(adonis_pvalue, file=paste("beta_pcoa_stat", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
    }
  }

  # loop for each group pair
  dis_table = as.matrix(dis_mat)
  if (pairwise) {
    compare_data = as.vector(unique(metadata$group))
    len_compare_data = length(compare_data)
    for(i in 1:(len_compare_data-1)) {
      for(j in (i+1):len_compare_data) {
        tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
        print(tmp_compare)
        da_adonis(tmp_compare)
      }
    }
  }else {
    compare_data = read.table(pairwise_list, sep="\t", check.names=F, quote='', comment.char="")
    colnames(compare_data) = c("sampA", "sampB")
    for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
  }
  adonis_pvalue = date()
  write.table(adonis_pvalue, file=paste("beta_pcoa_stat", ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)

}
