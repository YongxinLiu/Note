#' Dimensionality reduction through PCA
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param color A condition to color data points by e.g. "AGE"
#' @param shape A condition to shape data points by e.g. "SEX"
#' @param cx Component on the x-axis e.g. 1
#' @param cy Component on the y-axis e.g. 2
#' @param cz Component on the z-axis e.g. 3
#' @param n_neighbors Number of nearest neighbors
#' @param metric Distance function e.g. c("euclidean", "manhattan")
#' @param n_epochs Number of iterations
#' @param init Initial embedding using eigenvector e.g c("spectral", "random")
#' @param min_dist Determines how close points appear in the final layout
#' @param datatype Datatype to use e.g. c("logcpm", "relabu", "counts")
#' @return A list with a plotly object and summary table
#'
#' @examples
#' data_dir = system.file("extdata/MAE.rds", package = "animalcules")
#' toy_data <- readRDS(data_dir)
#' result <- dimred_umap(toy_data,
#'                       tax_level="genus",
#'                       color="AGE",
#'                       shape="DISEASE",
#'                       cx=1,
#'                       cy=2,
#'                       datatype="logcpm")
#' result$plot
#' 
#' @import umap
#' @import dplyr
#' @import scales
#' @import plotly
#' @import magrittr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
dimred_umap <- function(MAE,
                        tax_level,
                        color,
                        shape=NULL,
                        cx=1,
                        cy=2,
                        cz=NULL,
                        n_neighbors=15,
                        metric=c("euclidean", "manhattan"),
                        n_epochs=200,
                        init=c("spectral", "random"),
                        min_dist=0.1,
                        datatype=c("logcpm", "relabu", "counts")) {

    # Default variables
    metric <- match.arg(metric)
    init <- match.arg(init)
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
    
    # Custom Parameters
    umap.custom <- umap.defaults
    umap.custom$n_neighbors <- n_neighbors
    umap.custom$n_components <- max(cx, cy, cz)
    umap.custom$metric <- metric
    umap.custom$n_epochs <- n_epochs
    umap.custom$init <- init
    umap.custom$min_dist <- min_dist
    
    # Run UMAP
    umap.data <- umap(df, config=umap.custom, method="naive")
    df.umap <- umap.data$layout 
    colnames(df.umap) <- paste("C", seq_len(ncol(df.umap)), sep="")
    
    # Merge in covariate information
    if (!is.null(shape)) {
        df.umap.m <- merge(df.umap, 
                           sam_table[, c(color, shape), drop=FALSE], by=0, all=TRUE)
        # When shape is required
        # Bypass duplicate colnames if color == shape
        shape <- colnames(df.umap.m)[ncol(df.umap.m)] 
        df.umap.m[[shape]] <- as.factor(df.umap.m[[shape]])
        
    } else {
        df.umap.m <- 
            merge(df.umap, sam_table[, color, drop=FALSE], by=0, all=TRUE)
        shape <- 'shape' # Referenced by plotly later
        df.umap.m[[shape]] <- 1 # Constant results in omitting shape
    }
    
    # Plotly | Scatterplot
    if (is.null(cz)) {
        
        # 2D Plot
        p <- plot_ly(df.umap.m,
                     x = as.formula(paste("~C", cx, sep = "")),
                     y = as.formula(paste("~C", cy, sep = "")),
                     mode = "markers",
                     color = as.formula(paste("~", color, sep = "")),
                     symbol = as.formula(paste("~", shape, sep = "")),
                     type = "scatter",
                     text = df.umap.m$Row.names,
                     marker = list(size = 10))
    } else {
        
        # 3D Plot
        p <- plot_ly(df.umap.m,
                     x = as.formula(paste("~C", cx, sep = "")),
                     y = as.formula(paste("~C", cy, sep = "")),
                     z = as.formula(paste("~C", cz, sep = "")),
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
                     text = df.umap.m$Row.names,
                     marker = list(size = 6))
    }
    
    p$p <- NULL # To suppress a shiny warning
    
    return(list(plot=p))
}
