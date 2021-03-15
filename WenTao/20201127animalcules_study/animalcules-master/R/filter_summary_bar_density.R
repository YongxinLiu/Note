#' Data visualization by bar plot / density plot
#'
#' @param MAE A multi-assay experiment object
#' @param samples_discard The list of samples to filter
#' @param filter_type Either 'By Microbes' or 'By Metadata'
#' @param sample_condition Which condition to check e.g. 'SEX'
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' result <- filter_summary_bar_density(toy_data,
#'                                  samples_discard = c('subject_2', 'subject_4'),
#'                                  filter_type = 'By Metadata',
#'                                  sample_condition = 'SEX')
#' result
#'
#' @import dplyr
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
filter_summary_bar_density <- function(MAE, 
                                       samples_discard=NULL, 
                                       filter_type, 
                                       sample_condition) {
    
    # Subset the data
    MAE_subset <- 
    mae_pick_samples(MAE = MAE, discard_samples = samples_discard)
    
    # Extract data
    microbe <- MAE_subset[["MicrobeGenetics"]]
    # host <- MAE_subset[['HostGenetics']]
    sam_table <- as.data.frame(colData(microbe))  # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[, rownames(sam_table)]  # organism x sample
    
    # Add count summary data to sample table
    sam_table[, "Reads"] = colSums(counts_table[, rownames(sam_table)])
    sam_table[, "Taxnum"] = apply(counts_table, 2, function(x) sum(x >= 1))
    
    # select filter type
    if (filter_type == "Microbes") {
        cov <- "Taxnum"
    } else {
        cov <- sample_condition
    }
    
    # Use density plot if the variable has more than 8 
    # unique values Use bar plot if
    # the variable has less than 8 unique values
    num_levels <- length(unique(unlist(sam_table[, cov])))
    
    if (num_levels > 8 & num_levels/nrow(sam_table) >= 0.3 
        & !is.character(unlist(sam_table[, 
        cov]))) {
        fit <- density(unlist(sam_table[, cov]))
        num.density <- plot_ly(x = fit$x, y = fit$y, 
                type = "scatter", mode = "lines", 
                fill = "tozeroy") %>% 
            layout(title = cov, xaxis=list(title=cov), yaxis=list(title="Density"))
        num.density$elementId <- NULL
        return(num.density)
    } else {
        cat.df = data.frame(table(sam_table[, cov]))
        cat.bar <- plot_ly(x = cat.df$Var1, y = cat.df$Freq, 
                type = "bar", showlegend = FALSE) %>% 
                layout(title = cov, xaxis = list(tickmode = "array", 
                        showticklabels = TRUE, 
                categoryorder = "trace"), yaxis = list(title = "Frequency"))
        cat.bar$elementId <- NULL
        return(cat.bar)
    }
}
