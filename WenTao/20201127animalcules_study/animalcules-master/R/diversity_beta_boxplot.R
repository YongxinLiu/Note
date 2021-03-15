#' Beta diversity boxplot
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param input_beta_method bray, jaccard
#' @param input_select_beta_condition Which condition to group samples
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' p <- diversity_beta_boxplot(toy_data,
#'                             tax_level = 'genus',
#'                             input_beta_method = 'bray',
#'                             input_select_beta_condition = 'DISEASE')
#' p
#'
#' @import dplyr
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#' @import ape
#' @import GUniFrac
#' @export

diversity_beta_boxplot <- function(MAE, 
                                    tax_level, 
                                    input_beta_method, 
                                    input_select_beta_condition) {
    
    # Extract data
    microbe <- MAE[["MicrobeGenetics"]]  #double bracket subsetting is easier
    # host <- MAE[['HostGenetics']]
    tax_table <- as.data.frame(rowData(microbe))  # organism x taxlev
    sam_table <- as.data.frame(colData(microbe))  # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[, rownames(sam_table)]  # organism x sample
    
    # Sum counts by taxon level and return counts
    counts_table %<>% # Sum counts by taxon level
    upsample_counts(tax_table, tax_level)
    
    # change tax table size
    tax_table <- tax_table[,1:which(colnames(tax_table) %in% tax_level)]
    
    # generate beta diversity
    if (input_beta_method %in% c("bray", "jaccard")){
        # Then use vegdist from vegan to generate a bray distance object:
        dist.mat <- vegan::vegdist(t(counts_table), method = input_beta_method)
        dist.mat <- as.matrix(dist.mat)        
    } else {
        # unifrac
        # factorize each column
        tax_table[sapply(tax_table, is.character)] <- lapply(tax_table[sapply(tax_table, is.character)], 
                                               as.factor)
        # create formula
        frm = as.formula(paste0("~", paste(colnames(tax_table), collapse ="/")))
        
        # create phylo object
        tr <- as.phylo(frm, data = tax_table)
        
        # add branch length
        tr <- suppressWarnings(compute.brlen(tr))
        
        # root phylo
        tr <- suppressWarnings(root(tr,1,resolve.root = TRUE))

        
        # count table
        ct_table <- as.data.frame(t(counts_table))
        ct_table[sapply(ct_table, is.numeric)] <- lapply(ct_table[sapply(ct_table, is.numeric)], 
                                               as.integer)
        
        unifracs <- suppressWarnings(GUniFrac(ct_table,tr)$unifracs)
        dw <- unifracs[, , "d_1"]		# Weighted UniFrac
        du <- unifracs[, , "d_UW"]		# Unweighted UniFrac	
        if (input_beta_method == 'unweighted unifrac'){
            dist.mat <- du
        } else {
            dist.mat <- dw
        }
    }
    
    # change condition name
    colnames(sam_table)[which(colnames(sam_table) == 
    input_select_beta_condition)] <- "condition"
    
    dist.within.a <- c()
    dist.within.b <- c()
    dist.between <- c()
    for (i in seq_len(nrow(dist.mat))) {
        for (j in seq_len(nrow(dist.mat))) {
            if (sam_table$condition[i] == 
                unique(sam_table$condition)[1] & sam_table$condition[j] == 
                unique(sam_table$condition)[1]) {
                dist.within.a <- c(dist.within.a, dist.mat[i, j])
            } else if (sam_table$condition[i] == 
                unique(sam_table$condition)[2] & sam_table$condition[j] == 
                unique(sam_table$condition)[2]) {
                dist.within.b <- c(dist.within.b, dist.mat[i, j])
            } else {
                dist.between <- c(dist.between, dist.mat[i, j])
            }
        }
    }
    y.axis <- list(title = paste(input_beta_method, "Distance", sep = " "))
    p <- plot_ly(y = ~dist.within.a, type = "box", 
            name = paste("Within", unique(sam_table$condition)[1])) %>% 
        add_trace(y = ~dist.within.b, 
            name = paste("Within", unique(sam_table$condition)[2])) %>% 
        add_trace(y = ~dist.between, 
            name = "Between 2 conditions") %>% layout(yaxis = y.axis)
    p$elementId <- NULL  # To suppress a shiny warning
    return(p)
}