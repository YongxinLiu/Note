#' Categorize continuous variables
#'
#' @param sam_table A sample x condition dataframe
#' @param sample_condition Continuous variable to categorize
#' @param new_label Column name for categorized variable
#' @param nbins Auto select ranges for n bins/categories
#' @param bin_breaks Manually select ranges for bins/categories
#' @param bin_labels Manually label bins/categories
#' @return A list with an updated sample table and before/after plots
#'
#' @examples
#' library(SummarizedExperiment)
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' microbe <- MultiAssayExperiment::experiments(toy_data)[[1]]
#' samples <- as.data.frame(colData(microbe))
#' result <- filter_categorize(samples,
#'                             sample_condition = 'AGE',
#'                             new_label='AGE_GROUP',
#'                             bin_breaks=c(0,55,75,100),
#'                             bin_labels=c('Young','Adult','Elderly'))
#' result$sam_table
#' result$plot.unbinned
#' result$plot.binned
#'
#' @import plotly
#' @import magrittr
#'
#' @export
filter_categorize <- function(sam_table, 
                            sample_condition, 
                            new_label, 
                            nbins = NULL, 
    bin_breaks = c(), bin_labels = c()) {
    
    # Numeric vector according to covariate of interest
    unbinned <- as.numeric(unlist(sam_table[, sample_condition]))
    
    # Plot unbinned condition
    fit <- density(unbinned)
    p.unbinned <- plot_ly(x = fit$x, y = fit$y, 
        type = "scatter", mode = "lines", 
        fill = "tozeroy") %>% 
        layout(title = sample_condition, margin = list(l = 0, 
        r = 0, t = 30, b = 30))
    p.unbinned$elementId <- NULL
    
    # If bins are provided use those
    if (length(bin_breaks) >= 2) {
        nlabs <- length(bin_breaks) - 1
    } else {
        bin_breaks <- nbins
        nlabs <- nbins
    }
    # Ensure provided labels are the correct length
    if (length(bin_labels) >= 1) {
        stopifnot(length(bin_labels) == nlabs)
    }
    
    # Categorized condition
    binned <- cut.default(unbinned, bin_breaks, bin_labels)
    
    # Add into sam table
    sam_table[, new_label] <- binned
    
    # Plot Binned Condition
    p.binned <- 
    plot_ly(y = binned, type = "histogram") %>% 
    layout(title = new_label, 
    xaxis = list(title = "Frequency"), yaxis = list(title = "Bins"), 
    margin = list(l = 80))
    p.binned$elementId <- NULL
    
    return(list(sam_table = sam_table, 
            plot.unbinned = p.unbinned, 
            plot.binned = p.binned))
}
