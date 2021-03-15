#' Get alpha diversity using gini
#'
#' @param x A list of counts
#' @return A single value
#'
#' @examples
#' gini_simpson(seq_len(10))
#'
#' @export


# x: Species count vector
gini_simpson <- function(x) {

    # Simpson index
    lambda <- simpson_index(x)

    # Gini-Simpson diversity
    1 - lambda

}

