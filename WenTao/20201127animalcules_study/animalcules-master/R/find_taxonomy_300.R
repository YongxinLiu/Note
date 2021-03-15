#' Find the taxonomy for maximum 300 tids
#'
#' @param tids Given taxonomy ids
#' @return taxondata Data with the taxonomy information
#' @import rentrez
#' @import XML

#' @examples
#' taxonLevels <- find_taxonomy_300(tids=1200)
#'
#' @export
find_taxonomy_300 <- function(tids) {
    # @importFrom httr set_config
    # set_config(httr::config(http_version = 0))
    if (is.null(tids)) {
        return(NULL)
    }
    na.vec <- c()
    for (i in seq_len(length(tids))){
        if(is.na(tids[i])){
            na.vec <- c(na.vec, i)
        }
    }
    r_fetch <- entrez_fetch(db = "taxonomy", id = tids, rettype = "xml")
    dat <- xmlToList(r_fetch)
    taxonLevels <- lapply(dat, function(x) x$LineageEx)
    if(!is.null(na.vec)){
        for(i in seq_len(length(na.vec))){
        taxonLevels <- append(taxonLevels, list(NA), na.vec[i]-1)
        }
    }
    return(taxonLevels)
}
