#' Output biom
#'
#' @param MAE A multi-assay experiment object
#' @param path_to_output The folder to output biom file
#' @return A message
#'
#' @import biomformat
#' @import MultiAssayExperiment
#' @importFrom Matrix Matrix
#' 
#' @export

write_to_biom <- function(MAE, path_to_output){
    # Extract data
    microbe <- MAE[["MicrobeGenetics"]]  #double bracket subsetting is easier
    
    tax_table <- as.data.frame(SummarizedExperiment::rowData(microbe))
    sam_table <- as.data.frame(SummarizedExperiment::colData(microbe))
    counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[,
                    rownames(sam_table)]  # organism x sample
    
    counts_table_dge <- as(Matrix::Matrix(as.matrix(counts_table)), "dgeMatrix")
    
    y = make_biom(counts_table, sam_table, tax_table)
    
    
    metadata <- list()
    tax_order_7 <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    
    for (i in 1:nrow(tax_table)){
        metadata[[i]] <- as.character(tax_table[i,])
        if (length(metadata[[i]]) == 7){
            for (j in 1:7){
                metadata[[i]][j] <- paste0(tax_order_7[j], metadata[[i]][j])
            }
        } else if (length(metadata[[i]]) == 6){
            for (j in 1:7){
                if (j < 7){
                    metadata[[i]][j] <- paste0(tax_order_7[j], metadata[[i]][j])
                } else{
                    metadata[[i]] <- c(metadata[[i]], tax_order_7[7])
                }
                
            }
        } else{
            warning("Please check why the number of taxonomy levels is less than 6.")
        }
    }
    
    data_dir = system.file("extdata/biom_example.rds", package = "animalcules")
    xx = readRDS(data_dir)
    
    y_rows_new <- NULL
    for (i in 1:nrow(tax_table)){
        mm <- xx$rows[1,]
        mm[1,1] <- rownames(tax_table)[i]
        mm$metadata$taxonomy[[1]] <- metadata[[i]]
        rownames(mm) <- as.character(i)
        rownames(mm$metadata) <- as.character(i)
        y_rows_new <- rbind(y_rows_new, mm)
    }
    
    y$rows <- y_rows_new
    
    
    y_cols_new <- NULL
    for (i in 1:nrow(sam_table)){
        mm <- xx$columns[1,]
        mm[1,1] <- rownames(sam_table)[i]
        mm$metadata <- sam_table[i,]
        rownames(mm) <- as.character(i)
        rownames(mm$metadata) <- as.character(i)
        y_cols_new <- rbind(y_cols_new, mm)
    }
    
    y$columns <- y_cols_new
    
    write_biom(y, path_to_output)
    return("Successfully output the biom file!")
}
