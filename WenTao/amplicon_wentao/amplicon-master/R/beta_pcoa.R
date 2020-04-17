# 绘制beta多样性PCoA图+置信椭圆 Beta PCoA + stat ellipse
#
# This is the function named 'beta_pcoa'
# which draw PCoA scatter plot with stat ellipse, and reture a ggplot2 object
#
#' @title Plotting beta diversity scatter plot
#' @description Input alpha index and metadata, and manual set alpha index and metadata column names.
#' ggplot2 show PCoA with color and stat ellipse.
#' @param dis_mat distance matrix, typical output of usearch -beta_div,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @param ellipse stat ellipse, T or F.
#' @param label sample name showing, T or F.
#' @param PCo principle coordinate used, default 12, alternative 13, or 23.
#' @details
#' By default, returns beta PCoA coordinate
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray_curtis, unifrac}
#' \item{other used indices: unifrac_binary, euclidean, manhatten}
#' }
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso beta_cpcoa
#' @examples
#' # Set 3 parameters: dis_mat, metadata and groupID
#' beta_pcoa(dis_mat = beta_bray_curtis, metadata, "genotype")
#' # Set 2 parameters: dis_mat, metadata, and groupID as default "genotype"
#' beta_pcoa(dis_mat = beta_bray_curtis, metadata)
#' # Set 6 parameters: dis_mat, metadata, and groupID as using "site",
#' beta_pcoa(dis_mat = beta_unifrac, metadata, groupID = "site", ellipse = F, label = T, PCo = 13)
#' @export

# load("../data/beta_unifrac.rda")
# beta_pcoa(dis_mat = beta_unifrac, metadata, groupID = "site", ellipse = F, label = T, PCo = 13)

beta_pcoa <- function(dis_mat, metadata, groupID = "genotype", ellipse = T, label = F, PCo = 12) {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # groupID = "genotype"
  # ellipse = T
  # label = F
  # PCo = 12

  # 交叉筛选
  idx = rownames(metadata) %in% rownames(dis_mat)
  metadata = metadata[idx,]
  dis_mat = dis_mat[rownames(metadata), rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  # colnames(sampFile)[1] = "group"

  # PCoA
  pcoa = cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig = pcoa$eig
  points = cbind(points, sampFile[rownames(points),])
  colnames(points) = c("x", "y", "z","group")

  # 按1、2轴绘图
  if (PCo == 12){
    p = ggplot(points, aes(x=x, y=y, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按1、3轴绘图
  if (PCo == 13){
    p = ggplot(points, aes(x=x, y=z, color=group))  +
      labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按2、3轴绘图
  if (PCo == 23){
    p = ggplot(points, aes(x=y, y=z, color=group))  +
      labs(x=paste("PCo 2 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCo 3 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  p = p + geom_point(alpha=.7, size=2) + theme_classic() + theme(text=element_text(family="sans", size=7))
  # 是否添加置信椭圆
  if (ellipse == T){
    p = p + stat_ellipse(level=0.68)
  }
  # 是否显示样本标签
  if (label == T){
    p = p + geom_text_repel(label = paste(rownames(points)), colour="black", size=3)
  }
  p
}
