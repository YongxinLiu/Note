







#--------两两差异分析-备用方法-注释掉

# beta_disper <- function(physeq,grouping_column,pvalue.cutoff=0.05,which_distance){
#   
#   meta_table <- data.frame(sample_data(physeq))
#   meta_table$Groups <- meta_table[,grouping_column]
#   # compute beta dispersion
#   mod<- vegan::betadisper(phyloseq::distance(physeq,method=which_distance),meta_table$Groups,type="centroid")
#   # compute pairwise beta dispersion for all levels in the grouping variable
#   pmod <- vegan::permutest(mod, permutations = 99, pairwise = TRUE)
#   p.values <- pmod$pairwise$observed
#   # extract significantly dispersed pairs and assign significant labels
#   p.values <- p.values[!is.na(p.values <= pvalue.cutoff)]
#   signi_label <- paste(cut(p.values,breaks=c(-Inf,0.001,0.01,0.05, Inf), label=c("***", "**", "*", ".")))
#   groups_compared <- names(p.values)
#   
#   betadisper_res <- data.frame(groups_compared, p.values, signi_label)
#   out <- list("betadisper_res"=betadisper_res, "pmod"=pmod)
#   return(out)
# }
# 
# 
# pairtest <- beta_disper(ps1_rela, "Group", pvalue.cutoff, dist)
# pairResult <- pairtest$betadisper_res

#---------总全部的群落差异分析-------备用方法，注释掉------------------------------
# map = as.data.frame(sample_data(ps1_rela))
# unif <- phyloseq::distance(ps1_rela , method=dist, type="samples")
# ado<- vegan::adonis(unif  ~ map$Group,permutations = 999)
# a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
# R2 <- paste("adonis:R ",a, sep = "")
# b = as.data.frame(ado$aov.tab[6])[1,1]
# p_v = paste("p: ",b, sep = "")
# title = paste(R2," ",p_v, sep = "")






























# 
# betadisper_pw <- beta_disper(physeq, grouping_column, pvalue.cutoff, 
#                              which_distance)
# betadisper_pw <- betadisper_pw$betadisper_res
# 
# 
# 
# 
# 
# library(microbiomeSeq)
# ord.res <- ordination(physeq, which_distance = "bray", method = ord_meths[i], grouping_column = "Depth", 
#                       pvalue.cutoff = 0.05)
# 
# # 提取stress
# sol <- ord.res$solution
# #提取adonis检验结果
# adn_res <- ord.res$adonis_res
# #提取两两检验结果
# betadisper_res <- ord.res$betadispersion
# # 提取分组
# groups <- ord.res$groups
# 
# # 提取解释度
# ord.res$solution$
# 
# # 提取作图坐标
# ord_res <- data.frame(x = sol$points[, 1], y = sol$points[, 2], Groups = groups)
# 
# head(ord_res )
# 
# 
# p <- microbiomeSeq::plot.ordination(ord.res,  pvalue.cutoff = 0.05, show.pvalues = TRUE)
# 
# print(p)
# 
# 
# ordination.res = ord.res
# method = "NMDS"
# pvalue.cutoff = 0.05
# show.pvalues = T
# N = 5
# 
# extra_marginspace = 0.35
#   # 提取stress
#   sol <- ordination.res$solution
#   #提取adonis检验结果
#   adn_res <- ordination.res$adonis_res
#   #提取两两检验结果
#   betadisper_res <- ordination.res$betadispersion
#   # 提取分组
#   groups <- ordination.res$groups
#   
#   veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, 
#                               npoints = 100) {
#     theta <- (0:npoints) * 2 * pi/npoints
#     Circle <- cbind(cos(theta), sin(theta))
#     t(center + scale * t(Circle %*% chol(cov)))
#   }
#   
#   # 提取作图坐标
#   ord_res <- data.frame(x = sol$points[, 1], y = sol$points[, 2], Groups = groups)
#   
#   # 对每组的坐标求取均值
#   ord_res.mean = aggregate(ord_res[, 1:2], list(group = ord_res$Groups), 
#                            mean)
#   plot.new()
#   
#   ord <- ordiellipse(sol, groups, display = "sites", 
#                      kind = "se", conf = 0.95, label = T)
#   dev.off()
#   
#   df_ell <- data.frame()
#   for (g in levels(ord_res$Groups)) {
#     if (g != "" && (g %in% names(ord)) && all(eigen(ord[[g]]$cov)$values > 
#                                               0)) {
#       df_ell <- rbind(df_ell, cbind(as.data.frame(with(ord_res[ord_res$Groups == 
#                                                                  g, ], veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, 
#                                                                                        ord[[g]]$scale))), Groups = g))
#     }
#   }
#   colnames(df_ell) <- c("x", "y", "Groups")
#   gg_color_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]
#   }
#   cols = gg_color_hue(length(unique(ord_res$Groups)))
#   p <- ggplot2::ggplot(data = ord_res, aes(x, y, colour = Groups))
#   p <- p + ggplot2::geom_point(alpha = 0.5, size = 2)
#   p
#   p <- p + ggplot2::theme_bw()
#   p <- p + ggplot2::annotate("text", x = ord_res.mean$x, 
#                              y = ord_res.mean$y, label = ord_res.mean$group, size = 6, 
#                              colour = cols, family = "Courier", fontface = "bold", 
#                              alpha = 0.8, vjust = 0.3)
#   p
#   # 添加置信椭圆
#   p <- p + ggplot2::geom_path(data = df_ell, aes(x = x, y = y), 
#                               size = 1, linetype = 1, alpha = 0.3)
#   
#   if (method == "NMDS") {
#     stress.value <- sol$stress
#     stress.label <- paste("STRESS=", round(stress.value, 
#                                            4))
#     p <- p + ggplot2::annotation_custom(grob = textGrob(label = stress.label, 
#                                                         hjust = 0, gp = gpar(cex = 1.5, fontsize = 8)), ymin = max(ord_res$y), 
#                                         ymax = max(ord_res$y), xmin = extra_marginspace + 
#                                           max(ord_res$x), xmax = extra_marginspace + max(ord_res$x))
#     p <- p + xlab("NMDS1") + ylab("NMDS2")
#   }
#   else if (method == "PCoA") {
#     p <- p + xlab(paste("Dim1 (", sprintf("%.4g", 
#                                           sol$eig[1]), "%)", sep = "")) + ylab(paste("Dim2 (", 
#                                                                                      sprintf("%.4g", sol$eig[2]), "%)", sep = ""))
#   }
#   gt <- NULL
#   
#   if (!is.null(adn_res)) {
#     adn_pvalue <- adn_res[[1]][["Pr(>F)"]][1]
#     adn_rsquared <- round(adn_res[[1]][["R2"]][1], 
#                           3)
#     signi_label <- paste(cut(adn_pvalue, breaks = c(-Inf, 
#                                                     0.001, 0.01, 0.05, Inf), label = c("***", "**", 
#                                                                                        "*", ".")))
#     adn_res_format <- bquote(atop(atop("PERMANOVA", 
#                                        R^2 == ~.(adn_rsquared)), atop("p-value=" ~ 
#                                                                         .(adn_pvalue) ~ .(signi_label), phantom())))
#     if (adn_pvalue <= pvalue.cutoff) {
#       gt <- plot_adonis_res(p, ord_res, adn_res_format, 
#                             extra_marginspace)
#     }
#   }
#   anova_table <- NULL
#   if (!is.null(betadisper_res)) {
#     anova_table <- plot_betadisper(betadisper_res, show.pvalues, 
#                                    pvalue.cutoff, N)
#     th <- sum(anova_table$heights)
#   }
#   out <- p
#   if (!is.null(gt) && !is.null(anova_table)) {
#     out <- gridExtra::grid.arrange(gt, anova_table, heights = unit.c(unit(1, 
#                                                                           "null"), th))
#   }
#   else if (is.null(gt) && !is.null(anova_table)) {
#     out <- gridExtra::grid.arrange(p, anova_table, heights = unit.c(unit(1, 
#                                                                          "null"), th))
#   }
#   return(out)
# }




