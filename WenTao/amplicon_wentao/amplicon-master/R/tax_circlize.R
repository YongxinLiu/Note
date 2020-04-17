# 样本或组的物种组成弦图 Circlize of taxonomy for samples and gorups
#
# This is the function named 'tax_circlize'
# which draw circle, and reture a circlize object
#
#' @title Plotting circlize of taxonomy for groups or samples
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
#' \item{group: group mean circlize}
#' \item{sample: each sample circlize}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso tax_circlize
#' @examples
#' # example data: feature table, rownames is OTU/taxonomy, colnames is SampleID
#' data(tax_phylum)
#' # example data: metadata or design, include SampleID, genotype and site
#' data(metadata)
#' # Set 4 parameters: set top 5 taxonomy, group by "genotype"
#' tax_circlize(tax_sum = tax_phylum, metadata, topN = 5, groupID = "genotype")
#' @export

load("../data/tax_phylum.rda")
load("../data/metadata.rda")

# tax_circlize(tax_sum = tax_phylum, metadata, topN = 5, groupID = "genotype")
tax_circlize <- function(tax_sum, metadata, topN = 5, groupID = "genotype") {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "reshape2", "circlize")
  for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
    suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  }

  # 测试默认参数
  # library(amplicon)
  # tax_sum = tax_phylum
  # topN = 5
  # groupID = "genotype"

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(tax_sum)
  metadata = metadata[idx,]
  tax_sum = tax_sum[, rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

#----按丰度降序排序#----
  mean_sort = as.data.frame(tax_sum[(order(-rowSums(tax_sum))), ])
  # 筛选前N类，其它归为Other，可设置不同组数
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort,other)
  rownames(mean_sort)[topN] = c("Other")
  # 保存变量备份，并输出至文件
  merge_tax=mean_sort

#----按组合并求均值#----

    # 转置样品名添加组名，并去除多余的两个样品列
    mat_t = t(merge_tax)
    mat_t2 = merge(sampFile, mat_t, by="row.names")
    mat_t2 = mat_t2[,c(-1)]

    # 按组求均值，转置，再添加列名
    mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
    df = do.call(rbind, mat_mean)[-1,]
    geno = mat_mean$group
    colnames(df) = geno

#----默认参数绘图-颜色随机#----

    # 设置图片文件名、长宽和字体大小
    pdf(file="circlize.pdf", width=89/25.4, height=89/25.4, pointsize=8)
    # 上方绘图和图例代码
    chordDiagram(df)
    # 绘图结束后写入文件
    dev.off()

#----指定颜色绘图+图例#----
    # 颜色设定
    library(RColorBrewer)
    grid.col = NULL
    # 分类学颜色，最多12种
    grid.col[rownames(df)] = brewer.pal(dim(df)[1], "Set1")
    # 定义分组颜色
    grid.col[colnames(df)] = brewer.pal(dim(df)[2], "Accent")

    pdf(file="circlize_legend.pdf", width=183/25.4, height=89/25.4, pointsize=8)
    # 上方绘图和图例代码
    chordDiagram(df, directional = TRUE,diffHeight = 0.03, grid.col = grid.col, transparency = 0.5)
    # 添加行-物种图例
    legend("left",pch=20,legend=rownames(df),col=grid.col[rownames(df)],bty="n",cex=1,pt.cex=3,border="black")
    # 添加列-分组图例
    legend("right",pch=20,legend=colnames(df),col=grid.col[colnames(df)],bty="n",cex=1,pt.cex=3,border="black")
    # 绘图结束后写入文件
    dev.off()
}

# 输出结果存在AI中编辑缺失字体问题？AdobePiStd
