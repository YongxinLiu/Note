# A notification ID
# id <- NULL

url <- a("website!", href="https://compbiomed.github.io/animalcules-docs/")
output$tab <- renderUI({
  tagList("Need help? Check the docs on our", url)
})

# reactive values shared thorough the shiny app
data_dir = system.file("extdata/MAE.rds", package = "animalcules")
vals <- reactiveValues(
    MAE = readRDS(data_dir),
    MAE_backup = MAE
)

observeEvent(input$upload_example,{
  withBusyIndicatorServer("upload_example", {
    if (input$example_data == "toy"){
      data_dir = system.file("extdata/MAE.rds", package = "animalcules")
    } else if (input$example_data == "tb"){
      data_dir = system.file("extdata/TB_example_dataset.rds", 
                             package = "animalcules")
    } else if (input$example_data == "asthma"){
      data_dir = system.file("extdata/asthma_example_dataset.rds", 
                             package = "animalcules")
    }
    MAE_tmp = readRDS(data_dir)
    vals$MAE <- MAE_tmp
    vals$MAE_backup <- MAE_tmp
    # Update ui
    update_inputs(session)
  })
})

update_inputs <- function(session) {
    updateCovariate(session)
    updateSample(session)
    updateTaxLevel(session)
    updateOrganisms(session)
}

updateCovariate <- function(session){
    MAE <- vals$MAE
    covariates <- colnames(colData(MAE))
    # use experience to choose categorical variables
    covariates.colorbar <- c()
    covariates_remove_index <- NULL
    for (i in seq_len(length(covariates))){
        num.levels <- length(unique(colData(MAE)[[covariates[i]]]))
        if (num.levels > 1 & num.levels < 6){
            if (num.levels/nrow(colData(MAE)) < 0.7){
                covariates.colorbar <- c(covariates.colorbar, covariates[i])
            }

        } else if(num.levels >= 6){
            if (num.levels/nrow(colData(MAE)) < 0.1){
                covariates.colorbar <- c(covariates.colorbar, covariates[i])
            }
            if (num.levels/nrow(colData(MAE)) == 1 & is.character(colData(MAE)[[covariates[i]]])){
              if (is.null(covariates_remove_index)){
                covariates_remove_index <- i
              }else{
                covariates_remove_index <- c(covariates_remove_index, i)
              }
            }
        }
    }
    if (!is.null(covariates_remove_index)){
      covariates <- covariates[-covariates_remove_index]
    }
    
    # choose the covariates that has 2 levels
    covariates.two.levels <- c()
    for (i in seq_len(length(covariates))){
        num.levels <- length(unique(colData(MAE)[[covariates[i]]]))
        if (num.levels == 2){
            covariates.two.levels <- c(covariates.two.levels, covariates[i])
        }
    }

    ## numeric cov
    sam_temp <- colData(MAE)
    num_select <- lapply(covariates, function(x) is_categorical(unlist(sam_temp[,x])))
    num_covariates <- covariates[!unlist(num_select)]

    # Filter
    updateSelectInput(session, "filter_type_metadata", choices = covariates)
    updateSelectInput(session, "filter_bin_cov", choices = num_covariates)

    # Relabu
    updateSelectInput(session, "relabu_bar_sample_conditions", choices = covariates)
    updateSelectInput(session, "relabu_bar_group_conditions", choices = c("ALL", covariates))
    updateSelectInput(session, "relabu_heatmap_conditions", choices = covariates)
    updateSelectInput(session, "relabu_box_condition", choices = covariates.colorbar)

    # Dimred
    updateSelectInput(session, "dimred_pca_color", choices = covariates)
    updateSelectInput(session, "dimred_pca_shape", choices = c("None", covariates.colorbar))
    updateSelectInput(session, "dimred_pcoa_color", choices = covariates)
    updateSelectInput(session, "dimred_pcoa_shape", choices = c("None", covariates.colorbar))
    updateSelectInput(session, "dimred_umap_color", choices = covariates)
    updateSelectInput(session, "dimred_umap_shape", choices = c("None", covariates.colorbar))
    updateSelectInput(session, "dimred_tsne_color", choices = covariates)
    updateSelectInput(session, "dimred_tsne_shape", choices = c("None", covariates.colorbar))

    # Diversity
    updateSelectInput(session, "select_alpha_div_condition", choices = covariates)
    updateSelectInput(session, "select_beta_condition", choices = covariates.two.levels)
    updateSelectInput(session, "bdhm_select_conditions", choices = covariates.colorbar)

    # Differential
    updateSelectInput(session, "da_condition", choices = covariates)
    updateSelectInput(session, "da_condition_covariate", choices = covariates)
    updateSelectInput(session, "da_condition_limma", choices = covariates)
    updateSelectInput(session, "da_condition_covariate_limma", choices = covariates)
    # Biomarker
    updateSelectInput(session, "select_target_condition_biomarker", choices = covariates.two.levels)

}

# update taxonomy levels
updateTaxLevel <- function(session){
    MAE <- vals$MAE
    tax.name <- colnames(rowData(MAE[['MicrobeGenetics']]))

    # Filter Dashboard
    updateSelectInput(session, "assay_count_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "assay_ra_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "assay_logcpm_taxlev", choices = tax.name, selected=tax.default)
    
    # Relabu
    updateSelectInput(session, "relabu_bar_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "relabu_heatmap_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "relabu_box_taxlev", choices = tax.name, selected=tax.default)

    # Diversity
    updateSelectInput(session, "taxl.alpha", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "taxl.beta", choices = tax.name, selected=tax.default)

    # Dim Reduction
    updateSelectInput(session, "dimred_pca_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "dimred_pcoa_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "dimred_umap_taxlev", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "dimred_tsne_taxlev", choices = tax.name, selected=tax.default)

    # Differential
    updateSelectInput(session, "taxl.da", choices = tax.name, selected=tax.default)
    updateSelectInput(session, "taxl.da_limma", choices = tax.name, selected=tax.default)
    

    # Biomarker
    updateSelectInput(session, "taxl_biomarker", choices = tax.name, selected="genus")
}

# update samples
updateSample <- function(session){
    MAE <- vals$MAE
    sam.name <- rownames(colData(MAE[['MicrobeGenetics']]))

    # Filter
    updateSelectInput(session, "filter_sample_dis", choices = sam.name)

    # Relabu
    updateSelectInput(session, "relabu_bar_sample_iso", choices = sam.name)
    updateSelectInput(session, "relabu_bar_sample_dis", choices = sam.name)
    updateSelectInput(session, "relabu_heatmap_sample_iso", choices = sam.name)
    updateSelectInput(session, "relabu_heatmap_sample_dis", choices = sam.name)
}

# update organisms
updateOrganisms <- function(session){
    MAE <- vals$MAE
    org.name <- rownames(as.data.frame(assays(MAE[['MicrobeGenetics']])))

    # Filter
    updateSelectInput(session, "filter_organism_dis", choices = org.name)
}

observeEvent(input$upload_animalcules,{
  withBusyIndicatorServer("upload_animalcules", {
    MAE_tmp <- readRDS(input$rdfile$datapath)
    #print(colData(MAE_tmp))
    vals$MAE <- MAE_tmp
    vals$MAE_backup <- MAE_tmp
    # Update ui
    update_inputs(session)
  })
})


## animalcules-id
observeEvent(input$upload_mae,{
  withBusyIndicatorServer("upload_mae", {
    MAE_list <- readRDS(input$rdfile_id$datapath)
    # check the length of the list
    if (length(MAE_list) == 1){
        MAE_tmp <- MAE_list[[1]]
    } else{
        if (input$mae_data_type == "em"){
        MAE_tmp <- MAE_list[['em']]
        } else{
        MAE_tmp <- MAE_list[['hit']]
        }
    }

    vals$MAE <- MAE_tmp
    vals$MAE_backup <- MAE_tmp
    # Update ui
    update_inputs(session)
  })
})

### BIOM
observeEvent(input$upload_biom,{
  withBusyIndicatorServer("upload_biom", {
  biom_obj <- biomformat::read_biom(input$biom_id$datapath)
  count_table <- as.data.frame(as(biomformat::biom_data(biom_obj), "matrix"))
  tax_table <- as.data.frame(as(biomformat::observation_metadata(biom_obj), "matrix"))
  metadata_table <- as.data.frame(as(biomformat::sample_metadata(biom_obj), "matrix"))
  
  # create MAE object
  se_mgx <-
      count_table %>%
      base::data.matrix() %>%
      S4Vectors::SimpleList() %>%
      magrittr::set_names("MGX")

  se_colData <-
      metadata_table %>%
      S4Vectors::DataFrame()

  se_rowData <-
      tax_table %>%
      base::data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      #dplyr::select(superkingdom, phylum, class, order, family, genus) %>%
      S4Vectors::DataFrame()

  microbe_se <-
      SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
      S4Vectors::SimpleList(MicrobeGenetics = microbe_se)

  MAE <-
      MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)
    vals$MAE <- MAE
    vals$MAE_backup <- MAE
    # Update ui
    update_inputs(session)
  })
})

observeEvent(input$uploadDataCount,{
  withBusyIndicatorServer("uploadDataCount", {

  count_table <- read.table(input$countsfile$datapath,
                       header = input$header.count,
                       row.names = 1,
                       stringsAsFactors = FALSE,
                       sep = input$sep.count,
                       comment.char="",
                       check.names = FALSE)

  tax_table <- read.table(input$taxon.table$datapath,
                            header = input$header.count,
                            sep = input$sep.count,
                            row.names= 1,
                            stringsAsFactors=FALSE,
                            comment.char="",
                            check.names = FALSE)

  metadata_table <- read.table(input$annotfile.count$datapath,
                            header = input$header.count,
                            sep = input$sep.count,
                            row.names=input$metadata_sample_name_col_count,
                            stringsAsFactors=FALSE,
                            comment.char="",
                            check.names = FALSE)


  # Choose only the samples in metadata that have counts data as well
    sample_overlap <- intersect(colnames(count_table), rownames(metadata_table))
    if (length(sample_overlap) < length(colnames(count_table))){
        print(paste("The following samples don't have metadata info:",
                    paste(colnames(count_table)[which(!colnames(count_table) %in% sample_overlap)],
                          collapse = ",")))
        count_table <- count_table[,which(colnames(count_table) %in% sample_overlap)]
    }
    metadata_table <- metadata_table[match(colnames(count_table), rownames(metadata_table)),,drop=FALSE]



  # Test and fix the constant/zero row
  row.remove.index <- c()
  if (sum(base::rowSums(as.matrix(count_table)) == 0) > 0){
      row.remove.index <- which(base::rowSums(as.matrix(count_table)) == 0)
      count_table <- count_table[-row.remove.index,]
  }
  
  species_overlap <- intersect(rownames(count_table), rownames(tax_table))
  if (length(species_overlap) < length(rownames(count_table))){
    print(paste("The following species don't have taxonomy info:",
    paste(rownames(count_table)[which(!rownames(count_table) %in% species_overlap)],
    collapse = ",")))
    count_table <- count_table[which(rownames(count_table) %in% species_overlap),]
  }
  tax_table <- tax_table[match(rownames(count_table), rownames(tax_table)), ]
  
  
  # replace spaces in tax name with underscore
  tax_table <- as.data.frame(apply(tax_table,
                                   2,
                                   function(x)gsub('\\s+', '_',x)))
  
  
  
  # create MAE object
  se_mgx <-
      count_table %>%
      base::data.matrix() %>%
      S4Vectors::SimpleList() %>%
      magrittr::set_names("MGX")

  se_colData <-
      metadata_table %>%
      S4Vectors::DataFrame()

  se_rowData <-
      tax_table %>%
      base::data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      #dplyr::select(superkingdom, phylum, class, order, family, genus) %>%
      S4Vectors::DataFrame()

  microbe_se <-
      SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
      S4Vectors::SimpleList(MicrobeGenetics = microbe_se)

  MAE <-
      MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)


  # update vals
  vals$MAE <- MAE
  vals$MAE_backup <- MAE


  # Update ui
  update_inputs(session)
  })
})

observeEvent(input$uploadDataCountTi,{
  withBusyIndicatorServer("uploadDataCountTi", {

  count_table <- read.table(input$countsfileTi$datapath,
                       header = input$header.countTi,
                       row.names = 1,
                       stringsAsFactors = FALSE,
                       sep = input$sep.countTi,
                       comment.char="",
                       check.names = FALSE)
  metadata_table <- read.table(input$annotfile.countTi$datapath,
                            header = input$header.countTi,
                            sep = input$sep.countTi,
                            row.names=input$metadata_sample_name_col_countTi,
                            stringsAsFactors=FALSE,
                            comment.char="",
                            check.names = FALSE)


  # Choose only the samples in metadata that have counts data as well
    sample_overlap <- intersect(colnames(count_table), rownames(metadata_table))
    if (length(sample_overlap) < length(colnames(count_table))){
        print(paste("The following samples don't have metadata info:",
                    paste(colnames(count_table)[which(!colnames(count_table) %in% sample_overlap)],
                          collapse = ",")))
        count_table <- count_table[,which(colnames(count_table) %in% sample_overlap)]
    }
    metadata_table <- metadata_table[match(colnames(count_table), rownames(metadata_table)),,drop=FALSE]



  # Test and fix the constant/zero row
  row.remove.index <- c()
  if (sum(base::rowSums(as.matrix(count_table)) == 0) > 0){
      row.remove.index <- which(base::rowSums(as.matrix(count_table)) == 0)
      count_table <- count_table[-row.remove.index,]
  }
  
  ids <- rownames(count_table)
  if (startsWith(ids[1],'ti|')){
      tids <- unlist(lapply(ids, FUN = grep_tid))
  }
  
  if (sum(is.na(tids)) > 0){
      tid_remove <- which(is.na(tids))
      ids <- ids[-tid_remove]
      tids <- tids[-tid_remove]
      count_table <- count_table[-tid_remove,]
  }

  taxonLevels <- find_taxonomy(tids)
  tax_table <- find_taxon_mat(ids, taxonLevels)


  # replace spaces in tax name with underscore
  tax_table <- as.data.frame(apply(tax_table,
                                  2,
                                  function(x)gsub('\\s+', '_',x)))
  
  species_overlap <- intersect(rownames(count_table), rownames(tax_table))
  if (length(species_overlap) < length(rownames(count_table))){
    print(paste("The following species don't have taxonomy info:",
    paste(rownames(count_table)[which(!rownames(count_table) %in% species_overlap)],
    collapse = ",")))
    count_table <- count_table[which(rownames(count_table) %in% species_overlap),]
  }
  tax_table <- tax_table[match(rownames(count_table), rownames(tax_table)), ]

 
  # create MAE object
  se_mgx <-
      count_table %>%
      base::data.matrix() %>%
      S4Vectors::SimpleList() %>%
      magrittr::set_names("MGX")

  se_colData <-
      metadata_table %>%
      S4Vectors::DataFrame()

  se_rowData <-
      tax_table %>%
      base::data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      dplyr::select(superkingdom, phylum, class, order, family, genus) %>%
      S4Vectors::DataFrame()

  microbe_se <-
      SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
      S4Vectors::SimpleList(MicrobeGenetics = microbe_se)

  MAE <-
      MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)

  # update vals
  vals$MAE <- MAE
  vals$MAE_backup <- MAE

  # Update ui
  update_inputs(session)
  })
})

observeEvent(input$uploadDataPs, {
  withBusyIndicatorServer("uploadDataPs", {


    df.path.vec <- c()
    df.name.vec <- c()
    for(i in seq_len(length(input$countsfile.pathoscope[,1]))){
        df.path.vec[i] <- input$countsfile.pathoscope[[i, 'datapath']]
        df.name.vec[i] <- input$countsfile.pathoscope[[i, 'name']]
    }

    datlist <- read_pathoscope_data(input_dir, pathoreport_file_suffix = input$report_suffix,
                                  use.input.files = TRUE,
                                  input.files.path.vec = df.path.vec,
                                  input.files.name.vec = df.name.vec)
    count_table <- datlist$countdat

    metadata_table <- read.csv(input$annotfile.ps$datapath,
                              header = input$header.ps,
                              sep = input$sep.ps,
                              row.names=input$metadata_sample_name_col,
                              stringsAsFactors=FALSE,
                              strip.white=TRUE)

    # Choose only the samples in metadata that have counts data as well
    sample_overlap <- intersect(colnames(count_table), rownames(metadata_table))
    if (length(sample_overlap) < length(colnames(count_table))){
        print(paste("The following samples don't have metadata info:",
                    paste(colnames(count_table)[which(!colnames(count_table) %in% sample_overlap)],
                          collapse = ",")))
        #shinyalert::shinyalert("Notice!", "Some samples miss", type = "info")

        count_table <- count_table[,which(colnames(count_table) %in% sample_overlap)]
    }
    metadata_table <- metadata_table[match(colnames(count_table), rownames(metadata_table)),,drop = FALSE]

    # print("read in done!")
    # Test and fix the constant/zero row
    row.remove.index <- c()
    if (sum(base::rowSums(as.matrix(count_table)) == 0) > 0){
        row.remove.index <- which(base::rowSums(as.matrix(count_table)) == 0)
        count_table <- count_table[-row.remove.index,]
    }

    ids <- rownames(count_table)
    tids <- unlist(lapply(ids, FUN = grep_tid))
    if (sum(is.na(tids)) > 0){
        tid_remove <- which(is.na(tids))
        ids <- ids[-tid_remove]
        tids <- tids[-tid_remove]
        count_table <- count_table[-tid_remove,]
    }

    taxonLevels <- find_taxonomy(tids)
    tax_table <- find_taxon_mat(ids, taxonLevels)

    
    # replace spaces in tax name with underscore
   tax_table <- as.data.frame(apply(tax_table,
                                    2,
                                    function(x)gsub('\\s+', '_',x)))
  
    
  # create MAE object
  se_mgx <-
      count_table %>%
      base::data.matrix() %>%
      S4Vectors::SimpleList() %>%
      magrittr::set_names("MGX")

  se_colData <-
      metadata_table %>%
      S4Vectors::DataFrame()

  se_rowData <-
      tax_table %>%
      base::data.frame() %>%
      dplyr::mutate_all(as.character) %>%
      dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
      S4Vectors::DataFrame()

  microbe_se <-
      SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
      S4Vectors::SimpleList(MicrobeGenetics = microbe_se)

  MAE <-
      MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)
  ### debug
  # saveRDS(MAE, "~/Desktop/test.RDS")

  # update vals
  vals$MAE <- MAE
  vals$MAE_backup <- MAE

  # Update ui
  update_inputs(session)


  })
})

# Data input summary
output$contents.count <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$countsfile.pathoscope)){
        if (input$uploadChoiceAdv == "pathofiles"){
        req(input$countsfile.pathoscope)
        df <- read.table(input$countsfile.pathoscope[[1, 'datapath']],
                       skip = 1,
                       header = TRUE,
                       sep = input$sep.ps,
                       comment.char="",
                       check.names = FALSE
                       )
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))


output$biom.count <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$biom_id)){
        if (input$uploadChoiceAdv == "biom"){
        req(input$biom_id)
        df <- as.data.frame(as(biom_data(read_biom(input$biom_id$datapath)), "matrix"))
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))


output$contents.meta <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$annotfile.ps)){
        if (input$uploadChoiceAdv == "pathofiles"){
        req(input$countsfile.pathoscope)

        df <- read.csv(input$annotfile.ps$datapath,
                       header = input$header.ps,
                       sep = input$sep.ps)
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))


output$biom.meta <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$biom_id)){
        if (input$uploadChoiceAdv == "biom"){
        req(input$biom_id)
        df <- as.data.frame(as(sample_metadata(read_biom(input$biom_id$datapath)), "matrix"))
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))


output$biom.tax <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$biom_id)){
        if (input$uploadChoiceAdv == "biom"){
        req(input$biom_id)
        df <- as.data.frame(as(observation_metadata(read_biom(input$biom_id$datapath)), "matrix"))
        return(df)
        }
    }
},
options = list(
    paging = TRUE, scrollX = TRUE, pageLength = 5
))


### data input summary
output$contents.count.2 <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$countsfile)){
    if (input$uploadChoice == "count"){
      req(input$countsfile)
      df <- read.csv(input$countsfile$datapath,
                     header = input$header.count,
                     sep = input$sep.count,
                     check.names = FALSE)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.meta.2 <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$annotfile.count)){
    if (input$uploadChoice == "count"){
      req(input$annotfile.count)
      df <- read.csv(input$annotfile.count$datapath,
                     header = input$header.count,
                     sep = input$sep.count,
                     check.names = FALSE)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

### data input summary
output$contents.count.2Ti <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$countsfileTi)){
    if (input$uploadChoiceAdv == "countTi"){
      req(input$countsfileTi)
      df <- read.csv(input$countsfileTi$datapath,
                     header = input$header.countTi,
                     sep = input$sep.countTi,
                     check.names = FALSE)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.meta.2Ti <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.

  if (!is.null(input$annotfile.countTi)){
    if (input$uploadChoiceAdv == "countTi"){
      req(input$annotfile.countTi)
      df <- read.csv(input$annotfile.countTi$datapath,
                     header = input$header.countTi,
                     sep = input$sep.countTi,
                     check.names = FALSE)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))

output$contents.taxonomy <- DT::renderDataTable({

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.


  if (!is.null(input$taxon.table)){
    if (input$uploadChoice == "count"){
      req(input$taxon.table)

      df <- read.csv(input$taxon.table$datapath,
                     header = input$header.count,
                     sep = input$sep.count)
      return(df)
    }
  }
},
options = list(
  paging = TRUE, scrollX = TRUE, pageLength = 5
))
