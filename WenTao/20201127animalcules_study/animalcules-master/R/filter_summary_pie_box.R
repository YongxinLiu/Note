#' Data visualization by pie chart / box plot
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
#' result <- filter_summary_pie_box(toy_data,
#'                                  samples_discard = c('subject_2', 'subject_4'),
#'                                  filter_type = 'By Microbes',
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
filter_summary_pie_box <- function(MAE, 
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
    # Use density plot if the variable has more than 
    # 8 unique values Use pie chart if
    # the variable has less than 8 unique values
    num_levels <- length(unique(unlist(sam_table[, cov])))
    if (num_levels > 8 & num_levels/nrow(sam_table) >= 
        0.3 & !is.character(unlist(sam_table[, 
        cov]))) {
        vec <- unlist(sam_table[, cov])
        hover.txt <- paste(rownames(sam_table), ", ", vec, sep="")
        num.scatter <- plotly::plot_ly(y = vec, jitter = 0.3, 
            pointpos = -1.8, boxpoints = "all", 
            hoverinfo = "text",
            text=hover.txt,
            marker = list(color = "rgb(7,40,89)"), 
            line = list(color = "rgb(7,40,89)"), 
            name = cov, type = "box") %>% 
            layout(title = cov, yaxis=list(title=cov))
        num.scatter$elementId <- NULL
        return(num.scatter)
    } else {
        cat.df = data.frame(table(sam_table[, cov]))
        cat.pie <- plotly::plot_ly(cat.df, labels = ~Var1, 
            values = ~Freq, type = "pie", 
            showlegend = FALSE) %>% layout(title = cov)
        cat.pie$elementId <- NULL
        return(cat.pie)
    }
}