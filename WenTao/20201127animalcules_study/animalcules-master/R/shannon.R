#' Get alpha diversity using shannon
#'
#' @param x A list of counts
#' @return A single value
#'
#' @examples
#' shannon(seq_len(10))
#'
#' @export

# x: Species count vector
shannon <- function(x) {

    # Ignore zeroes
    x <- x[x > 0]

    # Species richness (number of species)
    S <- length(x)

    # Relative abundances
    p <- x/sum(x)

    # Shannon index
    (-sum(p * log(p)))
}
