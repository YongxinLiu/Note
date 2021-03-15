#' Reads the data from PathoScope reports and returns a list of
#' final guess relative abundance and count data
#'
#' @param input_dir Directory where the tsv files from PathoScope are located
#' @param pathoreport_file_suffix PathoScope report files suffix
#' @param use.input.files whether input dir to pathoscope files
#' or directly pathoscope files
#' @param input.files.path.vec vector of pathoscope file paths
#' @param input.files.name.vec vector of pathoscope file names
#' @return List of final guess relative abundance and count data
#' @importFrom utils read.table
#' @export

read_pathoscope_data <- function(input_dir = ".", 
                                 pathoreport_file_suffix = "-sam-report.tsv", 
                                 use.input.files = FALSE, 
                                 input.files.path.vec = NULL, 
                                 input.files.name.vec = NULL) {
    if (use.input.files == FALSE) {
        if (input_dir == ".") {
            input_dir <- getwd()
        }
        pattern <- paste("*", pathoreport_file_suffix, sep = "")
        filenames <- list.files(input_dir, 
                                pattern = pattern, 
                                full.names = TRUE)
        input.files.name.vec <- list.files(input_dir)
    } else {
        filenames <- input.files.path.vec
    }
    ltbl <- lapply(filenames, read.table, skip = 1,
            header = TRUE, sep = "\t", nrow = 100,
            comment.char = "", check.names = FALSE)
    lgenomes <- lapply(ltbl, function(tbl) {
        return((tbl[, 1]))
        #return(levels(tbl[, 1]))
    })
    genomes <- unique(unlist(lgenomes))
    # genomes <- c(genomes, 'others')
    lfl <- lapply(filenames, readLines, n = 1)
    lnumReads <- unlist(lapply(lfl, function(fl) {
        return(as.numeric(strsplit(fl, "\t")[[1]][2]))
    }))
    samplenames <- unlist(lapply(input.files.name.vec, function(x) {
        return(strsplit(x, "-sam-report.tsv")[[1]])
    }))
    # print(samplenames)
    countdat <- matrix(0L, 
    nrow = length(genomes), ncol = length(samplenames))
    for (i in seq(length(samplenames))) {
        index.tmp <- match(ltbl[[i]][, 1], genomes)
        countdat[index.tmp, i] <- floor(ltbl[[i]][, 4])
    }
    # integer
    countdat <- data.frame(countdat)
    rownames(countdat) <- genomes
    colnames(countdat) <- samplenames
    # saveRDS(countdat, '~/Desktop/test.rds')
    return(list(countdat = countdat))
}