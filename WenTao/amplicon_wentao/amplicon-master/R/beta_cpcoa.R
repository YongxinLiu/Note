# 绘制beta多样性 Constrained PCoA图+置信椭圆 Beta Constrained PCoA + stat ellipse
#
# This is the function named 'beta_cpcoa'
# which draw Constrained PCoA scatter plot with stat ellipse, and reture a ggplot2 object
#
#' @title Plotting beta diversity scatter plot of Constrained PCoA
#' @description Input rarefied OTU table and metadata, and manual set beta distance type and metadata column names.
#' ggplot2 show CPCoA with color and stat ellipse.
#' @param otutab rarefied OTU table, typical output of usearch -otutab_norm or -otutab_rare,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param dis distance type for caucluate, default "bray", alternative "manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup , binomial, chao, cao or mahalanobis".
#' @param groupID column name for groupID.
#' @param ellipse stat ellipse, T or F.
#' @param label sample name showing, T or F.
#' @details
#' By default, returns beta CPCoA coordinate
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray, manhattan, euclidean, jaccard}
#' \item{other used indices: canberra, kulczynski, gower, altGower, morisita, horn, mountford, raup , binomial, chao, cao or mahalanobis}
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
#' # example data: OTU table, rownames is OTU_xxx, colnames is SampleID
#' data(otutab)
#' # example data: metadata or design, include SampleID, genotype and site
#' data(metadata)
#' # Set 2 parameters: otu table, metadata, and distance type and groupID using default "bray" and "genotype"
#' beta_cpcoa(otutab, metadata)
#' # Set 4 parameters: otu table, metadata, distance type as "jaccard" and groupID as "site"
#' beta_cpcoa(otutab, metadata, "jaccard", "site")
#' # You can found a two-dimension seperated groups, by watch out p is significantly?
#' # Set 6 parameters: otutab, metadata, and groupID as using "site",
#' beta_cpcoa(otutab, metadata, dis = "bray", groupID = "genotype", ellipse = T, label = T)
#' @export

# load("../data/metadata.rda")
# load("../data/otutab.rda")
# beta_cpcoa(otutab, metadata, dis = "bray", groupID = "genotype", ellipse = T, label = T)


beta_cpcoa <- function(otutab, metadata, dis = "bray", groupID = "genotype", ellipse = T, label = F) {

  # 依赖关系检测与安装
  p_list = c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # dis = "bray"
  # groupID = "genotype"
  # ellipse = T
  # label = F

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx,]
  otutab = otutab[, rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  if (length(unique(sampFile$group)) > 2){

  # 函数提取CCA中主要结果
  variability_table = function(cca){
    chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
    variability_table = cbind(chi, chi/chi[1])
    colnames(variability_table) = c("inertia", "proportion")
    rownames(variability_table) = c("total", "constrained", "unconstrained")
    return(variability_table)
  }

  # Constrained analysis OTU table by genotype
  capscale.gen = capscale(t(otutab) ~ group, data=sampFile, add=F, sqrt.dist=T, distance=dis)

  # ANOVA-like permutation analysis
  perm_anova.gen = anova.cca(capscale.gen, permutations = 1000, parallel = 4)

  # generate variability tables and calculate confidence intervals for the variance
  var_tbl.gen = variability_table(capscale.gen)
  eig = capscale.gen$CCA$eig
  variance = var_tbl.gen["constrained", "proportion"]
  p.val = perm_anova.gen[1, 4]

  # extract the weighted average (sample) scores
  points = as.data.frame(capscale.gen$CCA$wa)
  points = cbind(sampFile, points[rownames(points),])

  # plot CPCo 1 and 2
  p = ggplot(points, aes(x=CAP1, y=CAP2, color=group)) + geom_point(alpha=.7, size=2) +
    labs(x=paste("CCA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("CCA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID) +
    ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",format(p.val, digits=2),sep="")) +
    theme_classic() + theme(text=element_text(family="sans", size=7))
  # 是否添加置信椭圆
  if (ellipse == T){
    p = p + stat_ellipse(level=0.68)
  }
  # 是否显示样本标签
  if (label == T){
    p = p + geom_text_repel(label = paste(rownames(points)), colour="black", size=3)
  }
  p
  }else{
    print("Selected groupID column at least have 3 groups!!!")
  }
}
