#' Find the taxonomy for unlimited tids
#'
#' @param tids Given taxonomy ids
#' @return A list of taxon levels with information
#'
#' @examples
#' taxonLevels <- find_taxonomy(tids=1200)
#'
#' @import rentrez

#' @export
find_taxonomy <- function(tids) {
    # @importFrom httr set_config 
    # set_config(httr::config(http_version = 0))
    if (is.null(tids)) {
        return(NULL)
    }
    if (length(tids) <= 300){
        taxonLevels <- find_taxonomy_300(tids)
    } else{
        taxonLevels <- list()
        batch.num <- ceiling(length(tids)/300)
        for (i in seq_len(batch.num)){
            if (i == batch.num){
                tids.batch <- tids[((i-1)*300 + 1):length(tids)]
            }else{
                tids.batch <- tids[((i-1)*300 + 1):(i*300)]
            }
        taxonLevels <- c(taxonLevels, find_taxonomy_300(tids.batch))
        }
    }
    return(taxonLevels)
}
