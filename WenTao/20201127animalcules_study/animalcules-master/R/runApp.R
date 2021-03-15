#' Run animalcules shiny app
#'
#' @import assertthat
#' @import covr
#' @import lattice
#' @import DT
#' @importFrom shinyjs addClass
#' @return The shiny app will open
#'
#' @param dev Run the applicaiton in developer mode
#' 
#' @examples 
#' \dontrun{
#' run_animalcules()
#' }
#' @export
run_animalcules <- function(dev=FALSE) {
    appDir <- system.file("shiny", package="animalcules")
    if (appDir == "") {
        stop("Could not find myapp. Try re-installing `mypackage`.", 
            call. = FALSE)
    }
    if (dev) {
        options(shiny.autoreload=TRUE)
    }
    shiny::runApp(appDir, display.mode="normal")
}