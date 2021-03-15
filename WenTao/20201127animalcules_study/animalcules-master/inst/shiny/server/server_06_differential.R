
  run.deseq2 <- eventReactive(input$run_deseq2, {
    withBusyIndicatorServer("run_deseq2", {
    MAE = vals$MAE
    differential_abundance(MAE = MAE,
                          tax_level=input$taxl.da,
                          input_da_condition=input$da_condition,
                          input_da_condition_covariate=input$da_condition_covariate,
                          min_num_filter = input$da.count.cutoff,
                          input_da_padj_cutoff = input$da.padj.cutoff,
                          method = "DESeq2")
    })

  })


  output$DeSeq2Table.new <- DT::renderDataTable({
    run.deseq2()
  },
  options=dtopts)


  run.limma <- eventReactive(input$run_limma, {
    withBusyIndicatorServer("run_limma", {
    MAE = vals$MAE
    differential_abundance(MAE = MAE,
                          tax_level=input$taxl.da_limma,
                          input_da_condition=input$da_condition_limma,
                          input_da_condition_covariate=input$da_condition_covariate_limma,
                          min_num_filter = input$da.count.cutoff_limma,
                          input_da_padj_cutoff = input$da.padj.cutoff_limma,
                          method = "limma")
    })

  })


  output$limmaTable.new <- DT::renderDataTable({
    run.limma()
  },
  options=dtopts)

