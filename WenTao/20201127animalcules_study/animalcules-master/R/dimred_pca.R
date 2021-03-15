#' Dimensionality reduction through PCA
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param color A condition to color data points by e.g. "AGE"
#' @param shape A condition to shape data points by e.g. "SEX"
#' @param pcx Principal component on the x-axis e.g. 1
#' @param pcy Principal component on the y-axis e.g. 2
#' @param pcz Principal component on the z-axis e.g. 3
#' @param datatype Datatype to use e.g. c("logcpm", "relabu", "counts")
#' @return A list with a plotly object and summary table
#'
#' @examples
#' data_dir = system.file("extdata/MAE.rds", package = "animalcules")
#' toy_data <- readRDS(data_dir)
#' result <- dimred_pca(toy_data,
#'                      tax_level="genus",
#'                      color="AGE",
#'                      shape="DISEASE",
#'                      pcx=1,
#'                      pcy=2,
#'                      datatype="logcpm")
#' result$plot
#' result$table
#'
#' @import dplyr
#' @import scales
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
dimred_pca <- function(MAE,
                        tax_level,
                        color,
                        shape=NULL,
                        pcx=1,
                        pcy=2,
                        pcz=NULL,
                        datatype=c("logcpm", "relabu", "counts")) {

    # Default variables
    datatype <- match.arg(datatype)

    # Extract data
    microbe <- MAE[['MicrobeGenetics']]
    #host <- MultiAssayExperiment::experiments(MAE)[[2]]
    tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

    df <- counts_table %>%
        # Sum counts by taxon level
        upsample_counts(tax_table, tax_level) %>%
        # Choose data type
            {
            if (datatype == "relabu") {
                counts_to_relabu(.)
            } else if (datatype == "logcpm") {
                counts_to_logcpm(.)
            } else {         .
            }
            } %>%
            # Fix constant/zero row
            {
                if (sum(base::rowSums(as.matrix(.)) == 0) > 0){
                . <- .[-which(base::rowSums(as.matrix(.)) == 0),]
            } else {
                .
            }
    } %>%
        # Transpose
        t()

    # PCA
    df.prcomp <- stats::prcomp(df, scale=TRUE)
    # Principle Components
    df.pca <- df.prcomp$x
    # Importance
    df.imp <- t(summary(df.prcomp)$importance)

    # Merge in covariate information
    if (!is.null(shape)) {
        df.pca.m <- merge(df.pca, 
        sam_table[, c(color, shape), drop=FALSE], by=0, all=TRUE)
        # When shape is required
        # Bypass duplicate colnames if color == shape
        shape <- colnames(df.pca.m)[ncol(df.pca.m)] 
        df.pca.m[[shape]] <- as.factor(df.pca.m[[shape]])

    } else {
        df.pca.m <- 
        merge(df.pca, sam_table[, color, drop=FALSE], by=0, all=TRUE)
        shape <- 'shape' # Referenced by plotly later
        df.pca.m[[shape]] <- 1 # Constant results in omitting shape
    }

    # Plotly | Scatterplot
    if (is.null(pcz)) {

        # 2D Plot
        p <- plot_ly(df.pca.m,
                    x = as.formula(paste("~PC", pcx, sep = "")),
                    y = as.formula(paste("~PC", pcy, sep = "")),
                    mode = "markers",
                    color = as.formula(paste("~", color, sep = "")),
                    symbol = as.formula(paste("~", shape, sep = "")),
                    type = "scatter",
                    text = df.pca.m$Row.names,
                    marker = list(size = 10))
    } else {

        # 3D Plot
        p <- plot_ly(df.pca.m,
                    x = as.formula(paste("~PC", pcx, sep = "")),
                    y = as.formula(paste("~PC", pcy, sep = "")),
                    z = as.formula(paste("~PC", pcz, sep = "")),
                    mode = "markers",
                    color = as.formula(paste("~", color, sep = "")),
                    symbol = as.formula(paste("~", shape, sep = "")),
                    symbols = c("circle", 
                        "square", 
                        "diamond", 
                        "cross", 
                        "square-open", 
                        "circle-open", 
                        "diamond-open", 
                        "x"),
                    type = "scatter3d",
                    text = df.pca.m$Row.names,
                    marker = list(size = 6))
    }

    p$p <- NULL # To suppress a shiny warning

    # Formatting importance table
    colnames(df.imp) = c("Standard Deviation",
                         "Variance Explained",
                         "Cumulative Variance")

    df.imp[,"Standard Deviation"] <- signif(df.imp[,"Standard Deviation"], 3)

    # Show variance as a percentage
    df.imp[,2] <- scales::percent(as.numeric(df.imp[,2]))
    df.imp[,3] <- scales::percent(as.numeric(df.imp[,3]))

    # Reorder
    df.imp <- as.data.frame(df.imp)
    df.imp$PC <- rownames(df.imp)
    df.imp <- df.imp[,c(4,1,2,3)]

    return(list(plot=p, table=df.imp))
}
