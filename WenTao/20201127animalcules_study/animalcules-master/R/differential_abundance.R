#' Differential abundance analysis
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param input_da_condition Which condition is the target condition
#' @param input_da_condition_covariate Covariates added to linear function
#' @param min_num_filter Minimum number reads mapped to this microbe
#' @param input_da_padj_cutoff adjusted pValue cutoff
#' @param method choose between DESeq2 and limma
#' 
#' @return A output dataframe
#'
#' @examples
#' data_dir = system.file("extdata/MAE.rds", package = "animalcules")
#' toy_data <- readRDS(data_dir)
#' differential_abundance(toy_data,
#' tax_level="phylum",
#' input_da_condition=c("DISEASE"),
#' min_num_filter = 2,
#' input_da_padj_cutoff = 0.5,
#' method = "DESeq2")
#' 
#'
#' @importFrom limma lmFit eBayes topTable
#' @import DESeq2
#' @import MultiAssayExperiment
#'
#' @export
differential_abundance <- function(MAE,
                                    tax_level,
                                    input_da_condition = c(),
                                    input_da_condition_covariate = NULL,
                                    min_num_filter = 5,
                                    input_da_padj_cutoff = 0.05,
                                    method = "DESeq2") {


    ## tables from MAE
    microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
    # organism x taxlev
    tax_table <- 
        as.data.frame(SummarizedExperiment::rowData(microbe)) 
    # sample x condition
    sam_table <- 
        as.data.frame(SummarizedExperiment::colData(microbe)) 
    # organism x sample
    counts_table <- 
        as.data.frame(SummarizedExperiment::assays(microbe))[,
                        rownames(sam_table)] 

    if (method == "DESeq2"){
        # Sum counts by taxon level
        count_table_tax <- counts_table %>%
                            upsample_counts(tax_table, tax_level)
        colnames_tmp <- colnames(count_table_tax)
        count_table_tax <- t(apply(count_table_tax, 1, as.integer))
        colnames(count_table_tax) <- colnames_tmp
        # sam table
        sam_table %<>% df_char_to_factor()
        # build the deseq2 formula
        if (is.null(input_da_condition_covariate)){
            dds_formula = 
            stats::as.formula(paste("~",input_da_condition, sep = " "))
        } else{
            dds_formula = stats::as.formula(paste("~",
                                    paste(
                                    paste(input_da_condition_covariate,
                                    collapse = " + "),
                                    input_da_condition,
                                    sep = " + "),
                                    sep = " "))
        }
        # run DEseq2
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_table_tax,
                                                colData = sam_table,
                                                design = dds_formula)
        dds <- DESeq2::DESeq(dds)
    
    
        # filter microbes with less than min_num_filter
        keep <- base::rowSums(DESeq2::counts(dds)) >= min_num_filter
        dds <- dds[keep,]
        #print(resultsNames(dds))
    
        # check if the condition has multiple levels
        if (length(resultsNames(dds)) == 2 | 
            length(resultsNames(dds)) - 
            length(input_da_condition_covariate) == 2){
            res <- DESeq2::results(dds)
            # reorder the result
            res = res[base::order(res$padj, na.last=NA), ]
        
        
            # reformat for reporting
            if (nrow(res) != 0){
            sigtab = res[(res$padj < input_da_padj_cutoff), ]
            if (nrow(sigtab) == 0){
            as.matrix("No differentially abundant items found!")
            } else{
                sigtab = as(sigtab, "data.frame")
                sigtab$padj <- 
            as.numeric(formatC(sigtab$padj, format = "e", digits = 2))
                sigtab$pValue <- 
            as.numeric(formatC(sigtab$pvalue, format = "e", digits = 2))
                sigtab$log2FoldChange <- 
            as.numeric(formatC(sigtab$log2FoldChange, 
                        format = "e", digits = 2))
                sigtab$microbe <- rownames(sigtab)
                rownames(sigtab) <- seq_len(nrow(sigtab))
                sigtab %<>% select(microbe, padj, pValue, log2FoldChange)
        
        
                num.1 <- c()
                num.2 <- c()
                # transform label into 1 and 0
                label.vec.num = 
            as.character((sam_table %>% select(input_da_condition))[,1])
                label.vec.save <- unique(label.vec.num)
                label.vec.num[label.vec.num == unique(label.vec.num)[1]] <- 1
                label.vec.num[label.vec.num != 1] <- 0
                label.vec.num <- as.numeric(label.vec.num)
                for (i in seq_len(nrow(sigtab))){
            species.index <- 
            which(rownames(count_table_tax) == sigtab[i,1])
            num.1 <- c(num.1, 
            sum((count_table_tax[species.index,
                        which(label.vec.num == 1)] > 0)))
            num.2 <- 
            c(num.2, sum((count_table_tax[species.index,
                        which(label.vec.num == 0)] > 0)))
                }
        
                sigtab <- cbind(sigtab, num.1)
                sigtab <- cbind(sigtab, num.2)
        
        
                df.output.prevalence <- 
            percent(round((num.1 + num.2)/ncol(count_table_tax),4))
                sigtab <- cbind(sigtab, df.output.prevalence)
        
        
                colnames(sigtab)[ncol(sigtab)-2] <- label.vec.save[1]
                colnames(sigtab)[ncol(sigtab)-1] <- label.vec.save[2]
                colnames(sigtab)[ncol(sigtab)] <- "prevalence"
        
        
                foldChange <- c()
                for (i in seq_len(nrow(sigtab))){
                foldChange[i] <- 
            round((max(as.numeric(c((sigtab[i,6] / sum(label.vec.num == 0)),
                        (sigtab[i,5] / sum(label.vec.num == 1))))) /
                min(as.numeric(c((sigtab[i,6] / sum(label.vec.num == 0)),
                (sigtab[i,5] / sum(label.vec.num == 1)))))), digits = 2)
                }
                sigtab <- cbind(sigtab, foldChange)
                colnames(sigtab)[ncol(sigtab)] <- 
            "Group Size adjusted fold change"
        
        
                # total num
                num_total <- length(label.vec.num)
                sigtab[,5] <- paste0(sigtab[,5], "/", sum(label.vec.num == 1))
                sigtab[,6] <- paste0(sigtab[,6], "/", sum(label.vec.num == 0))
                
                # if y is numeric, make the output table easier
                if (is.numeric((sam_table %>% 
                        select(input_da_condition))[,1])){
                        sigtab <- sigtab[,c(1,2,3,4,7)]
                }
        
                return(sigtab)
        
            }
        
            }else{
                return(as.matrix("No differentially abundant items found!"))
            }
        }else{
            # for multiple levels, 
            # we need to combine results for each comparison
            sigtab <- NULL
            label.vec = as.character((sam_table %>% 
                        select(input_da_condition))[,1])
            combination_mat <- utils::combn(sort(unique(label.vec)), 2)
            #print(combination_mat)
            for (j in seq(ncol(combination_mat))){
            res <- DESeq2::results(dds, contrast=c(input_da_condition, 
                                                combination_mat[1,j], 
                                                combination_mat[2,j]))
            if (nrow(res) > 0){
                res = res[base::order(res$padj, na.last=NA), ]
                sigtab_tmp = res[(res$padj < input_da_padj_cutoff), ]
                if (nrow(sigtab_tmp) > 0){
                        sigtab_tmp = as(sigtab_tmp, "data.frame")
                        sigtab_tmp$padj <- 
                        as.numeric(formatC(sigtab_tmp$padj, 
                                    format = "e", digits = 2))
                        sigtab_tmp$pValue <- 
                        as.numeric(formatC(sigtab_tmp$pvalue, 
                                    format = "e", digits = 2))
                        sigtab_tmp$log2FoldChange <- 
                        as.numeric(formatC(sigtab_tmp$log2FoldChange, 
                                    format = "e", digits = 2))
                        sigtab_tmp$microbe <- rownames(sigtab_tmp)
                        rownames(sigtab_tmp) <- seq_len(nrow(sigtab_tmp))
                        sigtab_tmp %<>% select(microbe, 
                                                padj, 
                                                pValue, 
                                                log2FoldChange)
                        num.1 <- c()
                        num.2 <- c()
                        # transform label into 1 and 0
                        label.vec.num = as.character((sam_table %>% 
                                    select(input_da_condition))[,1])
            label.vec.num[label.vec.num == combination_mat[1,j]] <- 1
            label.vec.num[label.vec.num == combination_mat[2,j]] <- 0
            label.vec.num <- as.numeric(label.vec.num)
                        for (i in seq_len(nrow(sigtab_tmp))){
                        species.index <- 
                        which(rownames(count_table_tax) == sigtab_tmp[i,1])
                        num.1 <- c(num.1, 
                        sum((count_table_tax[species.index,
                        which(label.vec.num == 1)] > 0)))
                        num.2 <- c(num.2, 
                        sum((count_table_tax[species.index,
                        which(label.vec.num == 0)] > 0)))
                        }
                        sigtab_tmp <- cbind(sigtab_tmp, num.1)
                        sigtab_tmp <- cbind(sigtab_tmp, num.2)
                        df.output.prevalence <- 
                        percent(round((num.1 + num.2)/
                                        ncol(count_table_tax),4))
                        sigtab_tmp <- cbind(sigtab_tmp, df.output.prevalence)
                        colnames(sigtab_tmp)[ncol(sigtab_tmp)-2] <- 
                        "experiment"
                        colnames(sigtab_tmp)[ncol(sigtab_tmp)-1] <- "control"
                        colnames(sigtab_tmp)[ncol(sigtab_tmp)] <- "prevalence"
                    
                        foldChange <- c()
                        for (i in seq_len(nrow(sigtab_tmp))){
                        foldChange[i] <- 
                        round((max(as.numeric(c((sigtab_tmp[i,6] / 
                        sum(label.vec == combination_mat[2,j])),
                        (sigtab_tmp[i,5] / 
                        sum(label.vec == combination_mat[1,j]))))) /
                        min(as.numeric(c((sigtab_tmp[i,6] / 
                        sum(label.vec == combination_mat[2,j])),
                        (sigtab_tmp[i,5] / 
                        sum(label.vec == combination_mat[1,j])))))), 
                        digits = 2)
                        }
                        sigtab_tmp <- cbind(sigtab_tmp, foldChange)
                        colnames(sigtab_tmp)[ncol(sigtab_tmp)] <- 
                        "Group Size adjusted fold change"
                        # total num
                        num_total <- length(label.vec.num)
                        sigtab_tmp[,5] <- 
                        paste0(combination_mat[1,j], 
                        ": ", 
                        sigtab_tmp[,5], 
                        "/", 
                        sum(label.vec == combination_mat[1,j]))
                        sigtab_tmp[,6] <- 
                        paste0(combination_mat[2,j], 
                        ": ", 
                        sigtab_tmp[,6], 
                        "/", 
                        sum(label.vec == combination_mat[2,j]))
                        # group 
                        sigtab_tmp[,9] <- 
                        paste0(combination_mat[1,j], 
                        " vs. ", combination_mat[2,j])
                        colnames(sigtab_tmp)[9] <- "Contrast"
                        # combine
                        sigtab <- rbind(sigtab, sigtab_tmp)
                }
            }
            }
            return(sigtab)
            }
    }else if (method == "limma"){
            # Sum counts by taxon level
            count_table_tax <- counts_table %>%
                        upsample_counts(tax_table, tax_level) %>%
                        counts_to_logcpm()
        colnames_tmp <- colnames(count_table_tax)
        count_table_tax <- t(apply(count_table_tax, 1, as.integer))
        colnames(count_table_tax) <- colnames_tmp
        # sam table
        sam_table %<>% df_char_to_factor()
        # filter low count microbes
        count_table_tax <- 
        count_table_tax[base::rowSums(count_table_tax) >= log10(min_num_filter),]
        #print(rowSums(count_table_tax))
        #print(min_num_filter)
        #print(rowSums(count_table_tax) >= min_num_filter)
        #print(count_table_tax)
        if (is.null(input_da_condition_covariate)){
            dds_formula <- 
            stats::as.formula(paste("~",input_da_condition, sep = " "))
            design <- model.matrix(dds_formula, sam_table)
        } else{
            dds_formula <- stats::as.formula(paste("~",
                        paste(
                        paste(input_da_condition_covariate,
                        collapse = " + "),
                        input_da_condition,
                        sep = " + "),
                        sep = " "))
            design <- model.matrix(dds_formula, sam_table)
        }
        #print(design)
        #print(str(count_table_tax))
        fit <- limma::lmFit(count_table_tax, design)
        ebayes <- limma::eBayes(fit)
        sigtab <- limma::topTable(ebayes, 
                                    adjust = "BH",
                                    number = nrow(count_table_tax),
                                    p.value=input_da_padj_cutoff)
        if (nrow(sigtab) == 0){
            sigtab[1,1] <- "No differentially abundant items found!"
            colnames(sigtab) <- "result"
            return(sigtab)
        }
        colnames(sigtab)[which(colnames(sigtab) == "adj.P.Val")] <- "padj"
        colnames(sigtab)[which(colnames(sigtab) == "P.Value")] <- "pValue"
        sigtab <- sigtab[,which(colnames(sigtab) %in% c("padj", "pValue"))]
        sigtab$microbe <- rownames(sigtab)
        rownames(sigtab) <- seq_len(nrow(sigtab))
        sigtab %<>% select(microbe, padj, pValue)
        return(sigtab) 
    }
}
