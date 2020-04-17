

#-----------统计大量的alpha多样性
library("microbiome")
# 来自microbiomedealpha多样性指标--这样一共包含25种多样性
ps = psRe
library("microbiome")
alp_mic = alpha(ps,index='all')
dim(alp_mic)
head(alp_mic)
colnames(alp_mic)
#phyloseq包中的alpha多样性,这里面的microbiome都有，所以就不用
ap_phy = estimate_richness(ps, measures =c("ACE","se.ACE"))
head(ap_phy)
alpfin = cbind(ap_phy,alp_mic)
head(alpfin)
finalp = cbind(alpfin,Richness)










plot_alpha <- function (physeq, x = "samples", color = NULL, shape = NULL, 
                        title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL, index = "all",
                        sortby = NULL) 
{
  
  
  erDF1 = estimate_richness(physeq, split = TRUE, measures = measures)
  
  tab <- alpha(physeq, index = index)
  
  
  erDF = cbind(erDF1,tab)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i, 
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures, 
                 ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, 
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby, 
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE, 
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color, 
                            shape = shape)
  p = ggplot(mdf, richness_map) + geom_point(na.rm = TRUE) + 
    stat_compare_means(comparisons=my_comparisons,label = "p.signif")+
    stat_compare_means()
  if (any(!is.na(mdf[, "se"]))) {
    # p = p + geom_errorbar(aes(ymax = value + se, ymin = value - 
    #                             se), width = 0.1)
  }
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, 
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}