#' Upsample a counts table to a higher taxon level
#'
#' @param counts_table A organism x sample data frame of counts
#' @param tax_table A organism x taxlev data frame of labels
#' @param higher_level Higher taxon level to upsample to
#' @return A organism x sample data frame of counts
#'
#' @examples	
#' toy_data <- readRDS(system.file("extdata/toy_data.rds", package = "animalcules"))
#' tax_table <- toy_data$tax_table
#' sam_table <- toy_data$sam_table
#' counts_table <- toy_data$counts_table 
#' counts_table <- upsample_counts(counts_table, tax_table, "phylum")
#' 
#' @import magrittr
#' @import reshape2
#' @import SummarizedExperiment
#'
#' @export
upsample_counts <- function(counts_table, tax_table, higher_level) {
    counts_table$higher_level = tax_table[[higher_level]]
    counts_table <- reshape2::melt(counts_table, id.vars = "higher_level") %>%
    S4Vectors::aggregate(. ~ 
        variable + higher_level, ., sum) %>% 
        reshape2::dcast(higher_level ~ variable) %>% 
        as.data.frame()
    rownames(counts_table) <- counts_table$higher_level
    counts_table$higher_level <- NULL
    # remove others
    counts_table <- counts_table[which(rownames(counts_table) != "others"),]
    return(counts_table)
}

#' Covert a counts table to a relative abundances table
#'
#' @param counts_table A organism x sample data frame of counts
#' @return A organism x sample data frame of relative abundances
#'
#' @examples
#' counts_to_relabu(matrix(seq_len(12),4))
#'
#' @import magrittr
#' @import SummarizedExperiment
#'
#' @export
counts_to_relabu <- function(counts_table) {
    prop.table(as.matrix(counts_table), 2) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(counts_table)) %>% 
    magrittr::set_rownames(rownames(counts_table))
}

#' Covert a counts table to a relative abundances table
#'
#' @param counts_table A organism x sample data frame of counts
#' @return A organism x sample data frame of logcpm counts
#'
#' @examples
#' logcpm <- counts_to_logcpm(as.data.frame(matrix(seq_len(12),4)))
#'
#' @import magrittr
#' @import SummarizedExperiment
#'
#' @export
counts_to_logcpm <- function(counts_table) {
    vapply(as.data.frame(counts_table), 
           function(x) log10(x * 1e+06/sum(x) + 1),
           c(rep(1.0,nrow(counts_table)))) %>%
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(counts_table)) %>% 
    magrittr::set_rownames(rownames(counts_table))
}

#' Modify samples of multi-assay experiment object
#'
#' @param MAE A multi-assay experiment object
#' @param isolate_samples Isolate specific samples e.g. c('SAM_01', 'SAM_02')
#' @param discard_samples Discard specific samples e.g. c('SAM_01', 'SAM_02')
#' @return A multi-assay experiment object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' subset <- mae_pick_samples(toy_data, 
#' isolate_samples=c('subject_9', 
#' 'subject_14'))
#'
#' @import MultiAssayExperiment
#'
#' @export
mae_pick_samples <- function(MAE, 
                             isolate_samples = NULL, 
                             discard_samples = NULL) {
    # Isolate all of these samples
    if (!is.null(isolate_samples)) {
        MAE <- MAE[, isolate_samples, ]
    }
    # Discard all of these samples
    if (!is.null(discard_samples)) {
        id = rownames(colData(MAE))
        id_isolate = id[!id %in% discard_samples]
        MAE <- MAE[, id_isolate, ]
    }
    return(MAE)
}

#' Modify organisms of multi-assay experiment object
#'
#' @param MAE A multi-assay experiment object
#' @param isolate_organisms Isolate specific organisms e.g. ti|001, ti|002
#' @param discard_organisms Discard specific organisms e.g. ti|001, ti|002
#' @return A multi-assay experiment object
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' subset <- mae_pick_organisms(toy_data, 
#' isolate_organisms=c('ti|001', 'ti|002'))
#'
#' @import MultiAssayExperiment
#'
#' @export
mae_pick_organisms <- function(MAE, 
                            isolate_organisms = NULL, 
                            discard_organisms = NULL) {
    # Isolate all of these organisms
    if (!is.null(isolate_organisms)) {
        MAE <- MAE[isolate_organisms, , ]
    }
    # Discard all of these organisms
    if (!is.null(discard_organisms)) {
        microbe <- MAE[["MicrobeGenetics"]]
        id = rownames(as.data.frame(assays(microbe)))
        id_isolate = id[!id %in% discard_organisms]
        MAE <- MAE[id_isolate, , ]
    }
    return(MAE)
}

#' Factorize all categorical columns
#'
#' @param df A sample x condition data frame
#' @return A sample x condition data frame
#'
#' @examples
#' df_char_to_factor(matrix(seq_len(12)))
#'
#'
#' @export
df_char_to_factor <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (typeof(df[, i]) == "character" | length(unique(df[, i])) < 4) {
            df[, i] <- as.factor(df[, i])
        }
    }
    return(df)
}

#' Format decimals to percentages
#'
#' @param x An array of decimals
#' @param digits number of digits
#' @param format f
#' @return An array of formatted strings
#'
#' @examples
#' nums <- c(0.42, 0.15, 0.4, 0.563, 0.2)
#' percent(nums)
#'
#' @export
percent <- function(x, digits = 2, format = "f") {
    paste0(formatC(100 * x, format = format, digits = digits), "%")
}

#' Check if object is categorical
#'
#' @param v A single value
#' @return Boolean
#'
#' @examples
#' nums <- 2
#' is_categorical(nums)
#'
#' @export
is_categorical <- function(v) {
    if (is.integer(v) || is.numeric(v)) {
        if (length(unique(v)) > 3) {
            return(FALSE)
        } else {
            return(TRUE)
        }
        
    } else {
        return(TRUE)
    }
}

#' check if integer(0)
#'
#' @param x A single value
#' @return Boolean
#'
#' @examples
#' nums <- 2
#' is_integer0(nums)
#'
#' @export
is_integer0 <- function(x) {
    is.integer(x) && length(x) == 0L
}

#' check if integer(1)
#'
#' @param x A single value
#' @return Boolean
#'
#' @examples
#' nums <- 2
#' is_integer1(nums)
#'
#' @export
is_integer1 <- function(x) {
    is.integer(x) && length(x) == 1L
}

#' Converts decimal percentage to string with specified digits
#'
#' @param v A single value
#' @param digits number of digits
#' @return Boolean
#'
#' @examples
#' nums <- 0.23
#' pct2str(nums)
#'
#' @export
pct2str <- function(v, digits = 2) {
    sprintf(paste0("%.", digits, "f"), v * 100)
}
