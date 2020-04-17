# 样本或组的物种组成堆叠柱状图 Stackplot of taxonomy for samples and gorups
#
# This is the function named 'tax_stackplot'
# which draw Constrained PCoA scatter plot with stat ellipse, and reture a ggplot2 object
#
#' @title Plotting stackplot of taxonomy for groups or samples
#' @description Input taxonomy composition, and metadata (SampleID and groupID). Then select top N high abundance taxonomy and group other low abundance. When Select samples can draw sample composition by facet groups. If used group can show mean of each group. Finally, return a ggplot2 object.
#' @param tax_sum composition matrix, like OTU table and rowname is taxonomy, typical output of usearch -sintax_summary;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param topN Top N taxonomy to show, default 8, alternative 4, 6, 10 ...;
#' @param groupID column name for groupID;
#' @param style group or sample, default group
#' @param sorted Legend sorted type, default abundance, alternative alphabet
#' @details
#' By default, returns top 8 taxonomy and group mean stackplot
#' The available style include the following:
#' \itemize{
#' \item{group: group mean stackplot}
#' \item{sample: each sample stackplot and facet by group}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso tax_stackplot
#' @examples
#' # example data: OTU table, rownames is OTU_xxx, colnames is SampleID
#' data(tax_phylum)
#' # example data: metadata or design, include SampleID, genotype and site
#' data(metadata)
#' # Set 2 parameters: tax_sum table, metadata, default include top 8 taxonomy, groupID is genotype, show group mean and sorted by abundance
#' tax_stackplot(tax_sum = tax_phylum, metadata)
#' # Set 4 parameters: set top10 taxonomy, group by "site", and group composition
#' tax_stackplot(tax_sum = tax_phylum, metadata, topN = 10, groupID = "site", style = "group", sorted = "abundance")
#' # Set 4 parameters: set top7 taxonomy, group by "genotype", and sample composition sorted by alphabet
#' tax_stackplot(tax_sum = tax_phylum, metadata, topN = 10, groupID = "genotype", style = "sample", sorted = "alphabet")
#' @export


# tax_stackplot(tax_sum = tax_phylum, metadata)

tax_stackplot <- function(tax_sum, metadata, topN = 8, groupID = "genotype", style = "group", sorted = "abundance") {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "reshape2")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # tax_sum = tax_phylum
  # topN = 8
  # groupID = "genotype"
  # style = "group"

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(tax_sum)
  metadata = metadata[idx,]
  tax_sum = tax_sum[, rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  # 按丰度降序排序
  mean_sort = as.data.frame(tax_sum[(order(-rowSums(tax_sum))), ])
  # 筛选前N类，其它归为Other，可设置不同组数
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort,other)
  rownames(mean_sort)[topN] = c("Other")
  # 保存变量备份，并输出至文件
  merge_tax=mean_sort

  if (style == "sample"){
  # 添加分类学并设置排序方式，默认字母，abundancer按丰度
  mean_sort$Taxonomy = rownames(mean_sort)
  data_all = as.data.frame(melt(mean_sort, id.vars=c("Taxonomy")))
  if (sorted == "abundance"){
  data_all$Taxonomy  = factor(data_all$Taxonomy, levels=rownames(mean_sort))
  }
  # set group facet
  data_all = merge(data_all, sampFile, by.x="variable", by.y = "row.names")

  # 按group分面，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  # 关闭x轴刻度和标签

  p = ggplot(data_all, aes(x=variable, y = value, fill = Taxonomy )) +
    geom_bar(stat = "identity",position="fill", width=1)+
    scale_y_continuous(labels = scales::percent) +
    facet_grid( ~ group, scales = "free_x", switch = "x") +
    theme(strip.background = element_blank())+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
    xlab("Groups")+ylab("Percentage (%)")+
    theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))+
    theme(text=element_text(family="sans", size=7))
  p
  }else{
    # 按组合并求均值

    # 转置样品名添加组名，并去除多余的两个样品列
    mat_t = t(merge_tax)
    mat_t2 = merge(sampFile, mat_t, by="row.names")
    mat_t2 = mat_t2[,c(-1)]

    # 按组求均值，转置，再添加列名
    mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
    mat_mean_final = do.call(rbind, mat_mean)[-1,]
    geno = mat_mean$group
    colnames(mat_mean_final) = geno
    mean_sort=as.data.frame(mat_mean_final)

    # 添加分类学并设置排序方式，默认字母，abundancer按丰度
    mean_sort$Taxonomy = rownames(mean_sort)
    data_all = as.data.frame(melt(mean_sort, id.vars=c("Taxonomy")))
    if (sorted == "abundance"){
      data_all$Taxonomy  = factor(data_all$Taxonomy, levels=rownames(mean_sort))
    }

    p = ggplot(data_all, aes(x=variable, y = value, fill = Taxonomy )) +
      geom_bar(stat = "identity",position="fill", width=0.7)+
      scale_y_continuous(labels = scales::percent) +
      xlab("Groups")+ylab("Percentage (%)")+ theme_classic()+
      theme(text=element_text(family="sans", size=7))
    p
  }
}
