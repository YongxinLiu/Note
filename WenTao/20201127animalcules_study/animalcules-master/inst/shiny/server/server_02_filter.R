## Filter by metadata
output$filter_metadata_params <- renderUI({

    MAE <- vals$MAE
    microbe <- MAE[['MicrobeGenetics']]
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    covdat <- sam_table[,input$filter_type_metadata]

    if (!is_categorical(covdat)) {
        sliderInput("filter_metadata_inp", "Include", min = min(covdat), max = max(covdat), value = c(min(covdat), max(covdat)))
    } else {
        selectizeInput("filter_metadata_inp", "Include", choices=unique(covdat), selected=unique(covdat), multiple=TRUE)
    }
})
observeEvent(input$filter_metadata_btn,{
    withBusyIndicatorServer("filter_metadata_btn", {

        MAE <- vals$MAE
        microbe <- MAE[['MicrobeGenetics']]
        sam_table <- as.data.frame(colData(microbe)) # sample x condition
        cov <- input$filter_type_metadata
        covdat <- sam_table[,cov]

        if (!is_categorical(covdat)) {
            minval <- input$filter_metadata_inp[1]
            maxval <- input$filter_metadata_inp[2]
            sam_table <- sam_table[covdat >= minval & covdat <= maxval,,drop=FALSE]
        } else {
            include <- input$filter_metadata_inp
            sam_table <- sam_table[covdat %in% include,,drop=FALSE]
        }
        samples <- rownames(sam_table)
        vals$MAE <- mae_pick_samples(MAE = vals$MAE, isolate_samples = samples)
        update_inputs(session)
    })
})

## Filter by microbes
# Filter by average read number
observeEvent(input$filter_microbes_read_btn,{
    withBusyIndicatorServer("filter_microbes_read_btn", {

        MAE <- vals$MAE
        microbe <- MAE[['MicrobeGenetics']]
        sam_table <- as.data.frame(colData(microbe)) # sample x condition
        counts_table <- as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

        minval <- input$filter_microbes_read_inp
        row_means <- apply(counts_table, 1, mean)
        organisms <- names(row_means[row_means >= minval])
        vals$MAE <- mae_pick_organisms(MAE, isolate_organisms = organisms)
        update_inputs(session)
    })
})
# Filter by average relative abundance
observeEvent(input$filter_microbes_rela_btn,{
    withBusyIndicatorServer("filter_microbes_rela_btn", {

        MAE <- vals$MAE
        microbe <- MAE[['MicrobeGenetics']]
        sam_table <- as.data.frame(colData(microbe)) # sample x condition
        counts_table <- as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample
        relabu_table <- counts_to_relabu(counts_table)

        minval <- input$filter_microbes_rela_min_inp
        maxval <- input$filter_microbes_rela_max_inp
        row_means <- apply(relabu_table, 1, mean)
        organisms <- names(row_means[row_means >= minval & row_means <= maxval])
        vals$MAE <- mae_pick_organisms(MAE, isolate_organisms = organisms)
        update_inputs(session)
    })
})
# Filter by average prevalence
observeEvent(input$filter_microbes_prev_btn,{
    withBusyIndicatorServer("filter_microbes_prev_btn", {

        MAE <- vals$MAE
        microbe <- MAE[['MicrobeGenetics']]
        sam_table <- as.data.frame(colData(microbe)) # sample x condition
        counts_table <- as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

        minval <- input$filter_microbes_prev_min_inp
        maxval <- input$filter_microbes_prev_max_inp
        row_means <- apply(counts_table, 1, function(x) (sum(x >= 1)/ncol(counts_table)))
        organisms <- names(row_means[row_means >= minval & row_means <= maxval])
        vals$MAE <- mae_pick_organisms(MAE, isolate_organisms = organisms)
        update_inputs(session)
    })
})

## Discard Samples and Organisms
observeEvent(input$filter_discard_btn,{
    withBusyIndicatorServer("filter_discard_btn", {
        MAE.1 <- mae_pick_samples(MAE = vals$MAE, discard_samples = input$filter_sample_dis)
        MAE.2 <- mae_pick_organisms(MAE = MAE.1, discard_organisms = input$filter_organism_dis)
        vals$MAE <- MAE.2
        update_inputs(session)
    })
})

## Reset data
observeEvent(input$filter_reset_btn,{
    withBusyIndicatorServer("filter_reset_btn", {
        vals$MAE <- vals$MAE_backup
        update_inputs(session)
    })
})

## Plots
output$filter_summary_top_plot <- renderPlotly({
    p <- filter_summary_bar_density(MAE = vals$MAE,
                                    samples_discard = c(),
                                    filter_type = input$filter_type,
                                    sample_condition = input$filter_type_metadata)
    return(p)
})

output$filter_summary_bottom_plot <- renderPlotly({
    p <- filter_summary_pie_box(MAE = vals$MAE,
                                samples_discard = c(),
                                filter_type = input$filter_type,
                                sample_condition = input$filter_type_metadata)
    return(p)
})

## Table
output$filter_summary_table <- renderTable({
    MAE = vals$MAE
    microbe <- MAE[['MicrobeGenetics']]
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    counts_table <- as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample
    relabu_table <- counts_to_relabu(counts_table)

    dat <- list()
    dat['Number of Samples'] <- round(nrow(sam_table))
    dat['Number of Covariates'] <- round(ncol(sam_table))
    dat['Number of Organisms'] <- round(nrow(counts_table))
    dat['Sample Mean Counts'] <- mean(apply(counts_table, 2, sum))
    dat['Sample Median Counts'] <- median(apply(counts_table, 2, sum))
    dat['Organism Mean Counts'] <- mean(apply(counts_table, 1, sum))
    dat['Organism Median Counts'] <- median(apply(counts_table, 1, sum))
    df <- as.data.frame(unlist(dat))

    # Formatting
    df$temp <- rownames(df)
    colnames(df) <- c("", "Summary Statistics")
    df <- df[,c(2,1)]
    df[,2] <- as.character(round(df[,2]))
    return(df)
})



## download
output$download_rds <- downloadHandler(filename = function() {
  paste("animalcules_data_", Sys.Date(), ".rds", sep="")
}, content = function(file) {
  saveRDS(vals$MAE, file=file)
})

output$download_biom <- downloadHandler(filename = function() {
  paste("animalcules_data_", Sys.Date(), ".biom", sep="")
}, content = function(file) {
  write_to_biom(vals$MAE, path_to_output=file)
})


## Categorize
output$filter_nbins <- renderUI({
    MAE <- vals$MAE
    microbe <- MAE[['MicrobeGenetics']]
    sam_table <- as.data.frame(colData(microbe)) # sample x condition
    vals <- unlist(sam_table[,input$filter_bin_cov,drop=TRUE])
    sliderInput("filter_nbins", label="Number of Bins", min=2, max=length(unique(vals)), value=2, step=1)
})
output$filter_bin_to1 <- renderPrint({
    x <- sort(as.numeric(unlist(strsplit(input$filter_bin_breaks,","))))
    print(x)
})
output$filter_bin_to2 <- renderPrint({
    x <- unlist(strsplit(input$filter_bin_labels,","))
    print(x)
})

output$filter_unbin_plot <- renderPlotly({
    MAE <- vals$MAE
    microbe <- MAE[['MicrobeGenetics']]
    samples <- as.data.frame(colData(microbe))
    result <- filter_categorize(samples,
                                sample_condition = input$filter_bin_cov,
                                new_label = input$filter_new_covariate)

    return(result$plot.unbinned)
})

do_categorize <- eventReactive(input$filter_create_bins, {
    MAE <- vals$MAE
    microbe <- MAE[['MicrobeGenetics']]
    samples <- as.data.frame(colData(microbe))

    nbins <- input$filter_nbins
    n <- input$filter_nbins

    # Overide custom bins if specified
    bin_breaks = sort(as.numeric(unlist(strsplit(input$filter_bin_breaks,","))))
    if (length(bin_breaks) > 1) {
      nbins = bin_breaks
      n = length(bin_breaks)-1
    }

    # Add custom labels only if the correct amount is sepecified
    bin_labels = unlist(strsplit(input$filter_bin_labels,","))
    fx_labels <- NULL
    if (length(bin_labels) == n) {
      fx_labels <- bin_labels
    }

    result <- filter_categorize(samples,
                                sample_condition = input$filter_bin_cov,
                                new_label = input$filter_new_covariate,
                                nbins = nbins,
                                bin_breaks = bin_breaks,
                                bin_labels = fx_labels)

    # Modify sample table
    MAE@colData <- DataFrame(result$sam_table)

    for (i in seq_len(length(MAE))){
        colData(MAE[[i]]) <- DataFrame(result$sam_table)
    }

    vals$MAE <- MAE
    print(colData(vals$MAE))
    update_inputs(session)

    return(result$plot.binned)
})

# Reaction to button pressing
output$filter_bin_plot <- renderPlotly({
    p <- do_categorize()
    return(p)
})




# Assays

# Render count table
find_assay_count <- eventReactive(input$view_assay_count, {
  MAE <- vals$MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe)) # sample x condition
  counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)] # organism x sample
  count_table_tax <- counts_table %>%
                    upsample_counts(tax_table, input$assay_count_taxlev)
  count_table_tax
})
output$assay_table_count <- DT::renderDataTable(find_assay_count(),
                                         options=dtopts, 
                                         rownames=TRUE)
output$download_assay_count <- downloadHandler(filename = function() {
  paste0("Assay_count_", input$assay_count_taxlev, ".csv", sep = "")
}, content = function(file) {
  df.out <- find_assay_count()
  write.csv(df.out, file)
})


# Render ra table
find_assay_ra <- eventReactive(input$view_assay_ra, {
  MAE <- vals$MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe)) # sample x condition
  counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)] # organism x sample
  relabu_table <- counts_table %>%
                upsample_counts(tax_table, input$assay_ra_taxlev) %>%
                counts_to_relabu() %>%
                base::as.data.frame()
  relabu_table
})
output$assay_table_ra <- DT::renderDataTable(find_assay_ra(),
                                         options=dtopts, 
                                         rownames=TRUE)
output$download_assay_ra <- downloadHandler(filename = function() {
  paste0("Assay_ra_", input$assay_ra_taxlev, ".csv", sep = "")
}, content = function(file) {
  df.out <- find_assay_ra()
  write.csv(df.out, file)
})

# Render logcpm table
find_assay_logcpm <- eventReactive(input$view_assay_logcpm, {
  MAE <- vals$MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe)) # sample x condition
  counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)] # organism x sample
  logcpm_table <- counts_table %>%
                upsample_counts(tax_table, input$assay_logcpm_taxlev) %>%
                counts_to_logcpm() %>%
                base::as.data.frame()
  logcpm_table
})
output$assay_table_logcpm <- DT::renderDataTable(find_assay_logcpm(),
                                         options=dtopts, 
                                         rownames=TRUE)
output$download_assay_logcpm <- downloadHandler(filename = function() {
  paste0("Assay_logcpm_", input$assay_logcpm_taxlev, ".csv", sep = "")
}, content = function(file) {
  df.out <- find_assay_logcpm()
  write.csv(df.out, file)
})


# Render taxonomy table
find_assay_tax <- eventReactive(input$view_assay_tax, {
  MAE <- vals$MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
  tax_table
})
output$assay_table_tax <- DT::renderDataTable(find_assay_tax(),
                                         options=dtopts, 
                                         rownames=TRUE)
output$download_assay_tax <- downloadHandler(filename = function() {
  "Assay_taxonomy.csv"
}, content = function(file) {
  df.out <- find_assay_tax()
  write.csv(df.out, file)
})


# Render annotation table
find_assay_annot <- eventReactive(input$view_assay_annot, {
  MAE <- vals$MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe)) # sample x condition
  sam_table
})
output$assay_table_annot <- DT::renderDataTable(find_assay_annot(),
                                         options=dtopts, 
                                         rownames=TRUE)
output$download_assay_annot <- downloadHandler(filename = function() {
  "Assay_annotation.csv"
}, content = function(file) {
  df.out <- find_assay_annot()
  write.csv(df.out, file)
})




