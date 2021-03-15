#' Alpha diversity statistical test
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param condition Which condition to group samples
#' @param alpha_metric Which alpha diversity metric to use
#' @param alpha_stat Which stat test to use
#' @return A dataframe
#'
#' @examples
#' data_dir = system.file("extdata/MAE.rds", package = "animalcules")
#' toy_data <- readRDS(data_dir)
#' p <- do_alpha_div_test(toy_data,
#'                        tax_level = "genus",
#'                        condition = "DISEASE",
#'                        alpha_metric = "shannon",
#'                        alpha_stat = "Wilcoxon rank sum test")
#' p
#'
#' @import dplyr
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#' @export

do_alpha_div_test <- function(MAE,
                            tax_level,
                            condition,
                            alpha_metric = c("inverse_simpson", "gini_simpson",
                                "shannon", "fisher", "coverage", "unit"),
                            alpha_stat = c("Wilcoxon rank sum test",
                                "T-test",
                                "Kruskal-Wallis")){

    # Extract data
    microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
    #host <- MAE[['HostGenetics']]
    tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

    # Sum counts by taxon level and return counts
    counts_table %<>%
        # Sum counts by taxon level
        upsample_counts(tax_table, tax_level)
    # calculate alpha diversity
    sam_table$richness <- diversities(counts_table, index = alpha_metric)
    colnames(sam_table)[ncol(sam_table)] <- "richness"
    colnames(sam_table)[which(colnames(sam_table) == condition)] <- "condition"
    # do the statistical test
    return(alpha_div_test(sam_table = sam_table,
                            alpha_stat = alpha_stat))
}


