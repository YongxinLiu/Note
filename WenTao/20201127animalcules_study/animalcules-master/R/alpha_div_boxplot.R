#' Alpha diversity boxplot
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param condition Which condition to group samples
#' @param alpha_metric Which alpha diversity metric to use
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' p <- alpha_div_boxplot(toy_data,
#'                        tax_level = 'genus',
#'                        condition = 'DISEASE',
#'                        alpha_metric = 'shannon')
#' p
#'
#' @import dplyr
#' @import plotly
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot labs
#' @importFrom stats model.matrix quantile t.test wilcox.test
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @export

alpha_div_boxplot <- function(MAE, 
                        tax_level, 
                        condition, 
                        alpha_metric = c("inverse_simpson", 
    "gini_simpson", "shannon", "fisher", "coverage", "unit")) {
    
    # Extract data
    microbe <- MAE[["MicrobeGenetics"]]  #double bracket subsetting is easier
    # host <- MAE[['HostGenetics']] organism x taxlev
    tax_table <- as.data.frame(SummarizedExperiment::rowData(microbe))
    # sample x condition
    sam_table <- as.data.frame(SummarizedExperiment::colData(microbe))
    counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[, 
                    rownames(sam_table)]  # organism x sample
    
    # Sum counts by taxon level and return counts
    counts_table %<>% upsample_counts(tax_table, tax_level)
    
    # calculate alpha diversity
    sam_table$richness <- diversities(counts_table, index = alpha_metric)
    colnames(sam_table)[ncol(sam_table)] <- "richness"
    colnames(sam_table)[which(colnames(sam_table) == condition)] <- "condition"
    
    
    # check if categorical variable
    if (!is.character(sam_table$condition)) {
        # plot alpha diversity boxplot
            g <- ggplot2::ggplot(sam_table, ggplot2::aes(condition, 
                        richness, 
                        text = rownames(sam_table), 
            color = condition)) + ggplot2::geom_point() + 
            ggplot2::labs(title = paste("Alpha diversity between ", 
            condition, " (", alpha_metric, ")", sep = ""))
    } else {
        # plot alpha diversity boxplot
        g <- ggplot2::ggplot(sam_table, ggplot2::aes(condition, 
                                    richness, 
                                    text = rownames(sam_table), 
            color = condition)) + ggplot2::geom_point() + 
            ggplot2::geom_boxplot() + 
            ggplot2::labs(title = paste("Alpha diversity between ", 
                        condition, " (", 
                        alpha_metric, ")", sep = ""))
    }
    
    # change y title
    g <- g + labs(y = alpha_metric)
    
    g <- ggplotly(g, tooltip = "text")
    g$elementId <- NULL  # To suppress a shiny warning
    return(g)
}