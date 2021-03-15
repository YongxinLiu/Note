#' Get alpha diversity
#'
#' @param sam_table A dataframe with 2 cols, richness and condition
#' @param alpha_stat Wilcoxon rank sum test or T-test for the test
#' @return A dataframe
#' @importFrom methods as
#' @importFrom stats as.formula cor.test density kruskal.test
#' @examples
#' df_test <- data.frame(richness = seq_len(10),
#' condition = c(rep(1,5), rep(0,5)))
#' alpha_div_test(df_test,alpha_stat='Wilcoxon rank sum test')
#'
#' @export
alpha_div_test <- function(sam_table, alpha_stat) {
    
    if (!is.character(sam_table$condition)) {
        tmp <- cor.test(sam_table$richness, sam_table$condition)
        output <- c(tmp$method, tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("Method", "P-value")
        return(output.table)
        
    }
    sam_table$condition <- as.factor(sam_table$condition)
    
    if (length(unique(sam_table$condition)) == 2) {
        tmp <- wilcox.test(richness ~ condition, data = sam_table)
        output <- c(tmp$p.value)
        output.table <- data.frame(output)
        rownames(output.table) <- c("P-value")
        colnames(output.table) <- tmp$method
        
        tmp <- t.test(richness ~ condition, data = sam_table)
        output_t <- c(tmp$p.value)
        output_t_table <- data.frame(output_t)
        rownames(output_t_table) <- c("P-value")
        colnames(output_t_table) <- tmp$method
        output.table <- cbind(output.table, output_t_table)
        return(output.table)
        
    } else if (length(unique(sam_table$condition)) == 3) {
        
        result.list <- list()
        sam_table.list <- list()
        for (i in seq_len(length(unique(sam_table$condition)))) {
            sam_table.list[[i]] <- 
            sam_table[which(sam_table$condition != 
                        unique(sam_table$condition)[i]), ]
            result.list[[i]] <- wilcox.test(richness ~ condition, 
                                    data = sam_table.list[[i]])
        }
        output.table <- NULL
        group.name <- c()
        for (i in seq_len(length(result.list))) {
            output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
            output.table <- cbind(output.table, output.tmp)
            group.name[i] <- paste(unique(sam_table.list[[i]]$condition), 
                                    collapse = " and ")
        }
        rownames(output.table) <- c("Method", "P-value")
        colnames(output.table) <- group.name
        
        result.list <- list()
        sam_table.list <- list()
        for (i in seq_len(length(unique(sam_table$condition)))) {
            sam_table.list[[i]] <- 
            sam_table[which(sam_table$condition != 
                        unique(sam_table$condition)[i]), 
                ]
            result.list[[i]] <- t.test(richness ~ condition, 
                                    data = sam_table.list[[i]])
        }
        output.table.t <- NULL
        group.name <- c()
        for (i in seq_len(length(result.list))) {
            output.tmp <- c(result.list[[i]]$method, result.list[[i]]$p.value)
            output.table.t <- cbind(output.table.t, output.tmp)
            group.name[i] <- paste(unique(sam_table.list[[i]]$condition), 
                                    collapse = " and ")
        }
        rownames(output.table.t) <- c("Method", "P-value")
        colnames(output.table.t) <- group.name
        
        tmp <- kruskal.test(richness ~ condition, data = sam_table)
        Overall <- c(tmp$method, tmp$p.value)
        output.table.k <- data.frame(Overall)
        rownames(output.table.k) <- c("Method", "P-value")
        output.table <- cbind(output.table.k, output.table)
        output.table <- cbind(output.table, output.table.t)
        return(output.table)
    } else {
        tmp <- kruskal.test(richness ~ condition, data = sam_table)
        output <- c(tmp$method, tmp$p.value)
        output.table.k <- data.frame(output)
        rownames(output.table.k) <- c("Method", "P-value")
        return(output.table.k)
    }
}
