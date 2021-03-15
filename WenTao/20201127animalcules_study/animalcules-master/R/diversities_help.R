#' Get alpha diversity
#'
#' @param counts_table A dataframe with organism x sample
#' @param index one of inverse_simpson,gini_simpson,shannon,fisher,coverage,unit
#' @param zeroes A boolean for whether to ignore zero values
#' @return A list of alpha diversity
#'
#' @examples
#' diversities_help(matrix(seq_len(12), nrow = 3),index='shannon')
#'
#' @export

diversities_help <- function(counts_table, index = "all", zeroes = TRUE) {
    
    if (length(index) > 1) {
        tab <- NULL
        for (idx in index) {
            tab <- cbind(tab, diversities_help(x, index = idx, zeroes = TRUE))
        }
        colnames(tab) <- index
        return(as.data.frame(tab))
    }
    
    if (index == "inverse_simpson") {
        ev <- apply(counts_table, 2, function(x) {
            inverse_simpson(x)
        })
    } else if (index == "gini_simpson") {
        ev <- apply(counts_table, 2, function(x) {
            gini_simpson(x)
        })
    } else if (index == "shannon") {
        ev <- apply(counts_table, 2, function(x) {
            shannon(x)
        })
    } else if (index == "unit") {
        ev <- apply(counts_table, 2, function(x) {
            sum(x>0)
        })        
    } else if (index == "fisher") {
        if (length(setdiff(unique(as.vector(counts_table)%%1), 0)) == 0) {
            ev <- vegan::fisher.alpha(counts_table, MARGIN = 2)
        } else {
            warning("Fisher diversity defined only for integers;
                the counts_table table contains non-integers. 
                Fisher not estimated.")
            ev <- NULL
        }
    } else if (index == "coverage") {
        ev <- unname(coverage(counts_table))
    }
    
    names(ev) <- colnames(counts_table)
    
    ev
    
}
