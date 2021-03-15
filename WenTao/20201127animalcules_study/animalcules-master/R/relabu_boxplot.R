#' Plot boxplots comparing different organism prevalence across conditions
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param condition Compare groups by condition e.g. 'SEX'
#' @param organisms Include organisms for plotting.
#' @param datatype counts, relative abundance,logcpm
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' p <- relabu_boxplot(toy_data,
#'                     tax_level='genus',
#'                     organisms=c('Escherichia', 'Actinomyces'),
#'                     condition='SEX',
#'                     datatype='logcpm')
#' p
#'
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
relabu_boxplot <- function(MAE, 
                            tax_level, 
                            condition, 
                            organisms = c(), 
                            datatype = c("counts", 
                                        "relative abundance", 
                                        "logcpm")) {
    
    # Default variables
    datatype <- match.arg(datatype)
    
    # Extract data
    microbe <- MAE[["MicrobeGenetics"]]
    # host <- MultiAssayExperiment::experiments(MAE)[[2]]
    tax_table <- as.data.frame(rowData(microbe))  # organism x taxlev
    sam_table <- as.data.frame(colData(microbe))  # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[, rownames(sam_table)]  # organism x sample
    
    # Ensure conditions are all factored
    sam_table %<>% df_char_to_factor()
    
    # Sum counts by taxon level and return counts
    df <- counts_table %>% # Sum counts by taxon level
    upsample_counts(tax_table, tax_level) %>% # Choose data type
    {
        if (datatype == "relative abundance") {
            counts_to_relabu(.)
        } else if (datatype == "logcpm") {
            counts_to_logcpm(.)
        } else {
            .
        }
    } %>% # Subset on chosen organisms
    .[organisms, , drop = FALSE] %>% # Transpose
    t() %>% # Back to data frame
    as.data.frame() %>% # Merge in covariate information
    merge(sam_table[, condition, drop = FALSE], 
        by = 0, all = TRUE) %>% # Melt for plotly
    reshape2::melt(by = organisms, variable.name = "organisms")
    
    p <- plot_ly(df, x = ~organisms, y = ~value, 
        color = ~df[, condition], type = "box") %>% 
        layout(boxmode = "group", xaxis = list(title = ""), 
            yaxis = list(title = datatype))
    p$p <- NULL  # To suppress a shiny warning
    return(p)
}
