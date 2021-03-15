#' Plot bar plots of sample and group level relative abundance
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param order_organisms A character list of organisms to send to top 
#' @param sort_by Sort bars by one of c("nosort", "conditions", "organisms", "alphabetically")
#' @param group_samples A bool specifying whether to group samples
#' @param group_conditions Group by one or more conditions e.g. "ALL" or "SEX"
#' @param sample_conditions Plot associatied conditions with samples.
#' @param isolate_samples Isolate specific samples e.g. c("SAM_01", "SAM_02")
#' @param discard_samples Discard specific samples e.g. c("SAM_01", "SAM_02")
#' @param show_legend A bool specifying whether or not to show organisms legend
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file("extdata/MAE.rds", package = "animalcules")
#' toy_data <- readRDS(data_dir)
#' p <- relabu_barplot(toy_data,
#'                     tax_level="family",
#'                     order_organisms=c('Retroviridae'),
#'                     sort_by="organisms",
#'                     sample_conditions=c('SEX', 'AGE'),
#'                     show_legend=TRUE)
#' p
#'
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
relabu_barplot <- function(MAE,
                            tax_level,
                            order_organisms=c(),
                            sort_by=c("nosort", "conditions", "organisms", "alphabetically"),
                            group_samples=FALSE,
                            group_conditions="ALL",
                            sample_conditions=c(),
                            isolate_samples=c(),
                            discard_samples=c(),
                            show_legend=TRUE) {

    # Default variables
    sort_by <- match.arg(sort_by)

    # Subset samples
    MAE <- mae_pick_samples(MAE, isolate_samples, discard_samples)

    # Extract data
    microbe <- MAE[['MicrobeGenetics']]
    #host <- MultiAssayExperiment::experiments(MAE)[[2]]
    tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

    # Ensure conditions are all factored
    sam_table %<>% df_char_to_factor()

    # Ensure sam table has the same subset of samples and is in correct order
    stopifnot(colnames(counts_table) == rownames(sam_table))

    # Sum counts by taxon level and return relative abundance
    relabu_table <- counts_table %>%
                    upsample_counts(tax_table, tax_level) %>%
                    counts_to_relabu() %>%
                    base::t() %>%
                    base::as.data.frame()

    # If grouping samples
    if (group_samples & !is.null(group_conditions)) {
        if (group_conditions == 'ALL') {
            relabu_table$covariate <- rep('ALL', nrow(relabu_table))
        } else {
            relabu_table$covariate <- sam_table[[group_conditions]]
        }

        # Group relative abundance by the covariate 
        # (mean relative abundance across organisms)
        relabu_table <- relabu_table %>%
                        reshape2::melt(id.vars = 'covariate') %>%
                        S4Vectors::aggregate(. ~variable+covariate , 
                                            ., mean) %>%
                        reshape2::dcast(formula = covariate~variable) %>%
                        magrittr::set_rownames(.[['covariate']]) %>%
                        dplyr::select(-one_of(c("covariate")))

        # Sam table becomes the grouped conditions
        sam_table <- rownames(relabu_table) %>%
                    as.data.frame() %>%
                    magrittr::set_colnames(c(group_conditions)) %>%
                    magrittr::set_rownames(rownames(relabu_table))
    }

    # Put selected organisms first
    relabu_table <- relabu_table[,order(colSums(relabu_table)),drop=FALSE]
    if (!is.null(order_organisms)) {
        org_order <- c(setdiff(colnames(relabu_table), 
                    order_organisms), 
                    rev(order_organisms))
        relabu_table <- relabu_table[,org_order]
    }
    
    # Put organisms alphabetically requested
    if (sort_by == "alphabetically") {
        org_order <- sort(colnames(relabu_table), decreasing=TRUE)
        relabu_table <- relabu_table[,org_order]
    }   
    
    # Order samples by organisms if not by conditons
    if (sort_by == "organisms") {
        for (i in seq_len(ncol(relabu_table))) {
            relabu_table <- relabu_table[order(relabu_table[,i]),]
        }
    }
    
    # If any conditions are selected make a side bar
    if (!is.null(sample_conditions) || 
        (group_samples && group_conditions != "ALL")) {
        if (!group_samples) {
            sam_table <- sam_table[,sample_conditions,drop=FALSE]
        }

        # Order samples by conditions if not by organisms
        if (sort_by == "conditions") {
            for (i in ncol(sam_table):1) {
                sam_table <- sam_table[order(sam_table[[i]]),,drop=FALSE]
            }
            # Reorder stacked barplot
            relabu_table <- 
            relabu_table[order(match(rownames(relabu_table), 
                rownames(sam_table))),,drop=FALSE]
        } else {
            sam_table <- 
            sam_table[order(match(rownames(sam_table), 
                    rownames(relabu_table))),,drop=FALSE]
        }

        if (nrow(sam_table) > 1) {
            # Retain hover-text information before conditions are factorized
            hover.txt <- c()
            for (i in seq_len(ncol(sam_table))) {
                hover.txt <- cbind(hover.txt, as.character(sam_table[[i]]))
            }

            # Normalized matrix
            mat <- sam_table %>%
                data.matrix() %>%
                apply(2, function(x)(x-min(x))/(max(x)-min(x)))

            # Plotly | Heatmap
            hm <- plotly::plot_ly(x = colnames(mat),
                                y = rownames(mat),
                                z = mat,
                                type = "heatmap",
                                showscale=FALSE,
                                hoverinfo = "x+y+text",
                                text=hover.txt) %>%
                                layout(xaxis = list(title = "",
                                        tickangle = -45),
                                        yaxis = list(showticklabels = FALSE,
                                                type = 'category',
                                                ticks = ""))
        }

    }

    # Plotly | Stacked Bar Plots
    relabu_table$samples <- rownames(relabu_table)

    sbp <- plotly::plot_ly(relabu_table,
                    y = ~samples,
                    x = relabu_table[[colnames(relabu_table)[1]]],
                    type = 'bar',
                    textposition= 'outside',
                    orientation = 'h',
                    name = substr(colnames(relabu_table)[1], 1, 40)) %>%
                    layout(font = list(size = 10),
                        xaxis = list(title = 'Relative Abundance',
                        automargin = TRUE),
                        yaxis = list(title = '',
                                    type = 'category',
                                    tickmode = "array",
                                    tickvals = rownames(relabu_table),
                                    showticklabels = FALSE,
                                    categoryorder = 'trace',
                                    automargin = TRUE),
                        barmode = 'stack',
                        showlegend = show_legend)
    for (i in 2:(ncol(relabu_table)-1)) {
        sbp <- 
        add_trace(sbp, x=relabu_table[[colnames(relabu_table)[i]]], 
            name=substr(colnames(relabu_table)[i], 1, 40))
    }

    # Create a multiplot if any conditions are selected
    if (exists("hm") && nrow(sam_table) > 1) {
        hm_sbp <- subplot(hm, sbp, widths=c(0.1,  0.9))
        hm_sbp$elementId <- NULL # To suppress a shiny warning
        return(hm_sbp)
    } else {
        sbp$elementId <- NULL # To suppress a shiny warning
        return(sbp)
    }
}
