# 绘制alpha多样性按组稀疏曲线 Alpha rarefracation curve
# This function named 'alpha_rare_curve',
# which draw curve and standard error with alpha_rare and metadata, and return a ggplot2 object

#' @title Plotting rarefracation curve for each group with anova statistics
#' @description Input alpha rare and metadata, and manual set metadata column names.
#' ggplot2 show lineplot, and standard error
#' @param alpha_rare alpha rarefracation matrix, typical output of usearch -alpha_div_rare;
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param groupID column name for groupID.
#' @details
#' By default, returns grouped curve, sample parameter to be continued
#' @return ggplot2 object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso alpha_boxplot
#' @examples
#' # Set three parameters: alpha_rare, metadata, groupID
#' alpha_rare_curve(alpha_rare, metadata, "genotype")
#' # Set two parameters: alpha_div, metadata, groupID as default "genotype"
#' alpha_rare_curve(alpha_rare, metadata)
#' # Set three parameters: alpha_rare, metadata, and groupID using site
#' alpha_rare_curve(alpha_rare, metadata, "site")
#' @export

# load("../data/alpha_rare.rda")
# alpha_rare_curve(alpha_rare, metadata, "site")




alpha_rare_curve <- function(alpha_rare, metadata, groupID = "genotype") {

  # 依赖关系检测与安装
  p_list = c("ggplot2","reshape2")
  for(p in p_list){
    if (!requireNamespace(p)){
    install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

  # 测试默认参数
  # groupID = "genotype"

  # 交叉筛选
  idx = rownames(metadata) %in% colnames(alpha_rare)
  metadata = metadata[idx,]
  alpha_rare = alpha_rare[,rownames(metadata)]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"

  # 求各组均值
  # 默认步长为1，折线不平滑，改为4减少锯齿
  df =alpha_rare[(1:25)*4,]
  # 转置df表格与实验设计合并，并去除第一列样品名
  mat_t = merge(sampFile, t(df), by="row.names")[,-1]
  # 按第一列合并求均值
  mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean)
  # 修正行名
  mat_mean_final = do.call(rbind, mat_mean)[-1,]
  geno = mat_mean$group
  colnames(mat_mean_final) = geno

  df=as.data.frame(round(mat_mean_final))
  df$x = rownames(df)
  df_melt = melt(df, id.vars=c("x"))

  # 求各组标准误
  # 转置df表格与实验设计合并，并去除第一列样品名
  se = function(x) sd(x)/sqrt(length(x)) # function for Standard Error
  mat_se = aggregate(mat_t[,-1], by=mat_t[1], FUN=se)
  mat_se_final = do.call(rbind, mat_se)[-1,]
  colnames(mat_se_final) = geno

  df_se=as.data.frame(round(mat_se_final))
  df_se$x = rownames(df_se)
  df_se_melt = melt(df_se, id.vars=c("x"))

  # 添加标准误到均值中se列
  df_melt$se=df_se_melt$value

  # 添加levels顺序，否则按字母顺序
  df_melt$x = factor(df_melt$x, levels=c(1:100))

  p = ggplot(df_melt, aes(x = x, y = value, group = variable, color = variable )) +
    geom_line()+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.5) +
    labs(x="Percentage (%)", y=paste("Richness"), color=groupID)+theme_classic()+
    scale_x_discrete(breaks = c(1:10)*10, labels = c(1:10)*10)+
    theme(text=element_text(family="sans", size=7))
  # xlab("Percentage (%)")+ylab("Richness")
  p
}
