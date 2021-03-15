#' Beta diversity test (by default we use bray-curtis distance)
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param input_beta_method bray, jaccard
#' @param input_select_beta_condition Which condition to group samples
#' @param input_select_beta_stat_method PERMANOVA,Kruskal-Wallis,Wilcoxon test
#' @param input_num_permutation_permanova number of permutations
#' @return A plotly object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' p <- diversity_beta_test(toy_data,
#'                          tax_level = 'genus',
#'                          input_beta_method = 'bray',
#'                          input_select_beta_condition = 'DISEASE',
#'                          input_select_beta_stat_method = 'PERMANOVA',
#'                          input_num_permutation_permanova = 999)
#' p
#'
#' @import dplyr
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#' @export

diversity_beta_test <- function(MAE, 
                                tax_level, 
                                input_beta_method, 
                                input_select_beta_condition, 
                                input_select_beta_stat_method, 
                                input_num_permutation_permanova = 999) {
    
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
        tr <- root(tr,1,resolve.root = TRUE)
        
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
        
        
        


    
    
    colnames(sam_table)[which(colnames(sam_table) == 
            input_select_beta_condition)] <- "condition"
    
    if (input_select_beta_stat_method == "PERMANOVA") {
        # bioconductor not allow set seed within R code set.seed(99)
        beta.div <- vegan::adonis2(dist.mat ~ condition, 
            data = sam_table, 
            permutations = input_num_permutation_permanova, 
            strata = "PLOT")
        beta.div
    } else {
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
                        unique(sam_table$condition)[2] & 
                    sam_table$condition[j] == 
                    unique(sam_table$condition)[2]) {
                    dist.within.b <- c(dist.within.b, dist.mat[i, j])
                } else {
                    dist.between <- c(dist.between, dist.mat[i, j])
                }
                
            }
        }
        dist.list <- list(dist.within.a, dist.within.b, dist.between)
        names(dist.list) <- 
        c(unique(sam_table$condition)[1], 
        unique(sam_table$condition)[2], 
        "between")
        if (input_select_beta_stat_method == "Wilcoxon rank sum test") {
            result.list <- list()
            group.name <- c()
            for (i in seq_len(length(dist.list))) {
                dist.list.tmp <- 
                dist.list[which(names(dist.list) != names(dist.list)[i])]
                group.name[i] <- paste(names(dist.list.tmp), collapse = " and ")
                result.list[[i]] <- 
                wilcox.test(dist.list.tmp[[1]], dist.list.tmp[[2]])
            }
            output.table <- NULL
            for (i in seq_len(length(result.list))) {
                output.tmp <- 
                c(result.list[[i]]$method, result.list[[i]]$p.value)
                output.table <- cbind(output.table, output.tmp)
            }
            rownames(output.table) <- c("Method", "P-value")
            colnames(output.table) <- group.name
            output.table
        } else {
            tmp <- 
            kruskal.test(list(dist.within.a, dist.within.b, dist.between))
            output <- c(tmp$method, tmp$p.value)
            output.table <- data.frame(output)
            rownames(output.table) <- c("Method", "P-value")
            output.table
        }
    }
}

