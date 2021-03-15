

  observeEvent(input$goButtonBiomarker, {
    withBusyIndicatorServer("goButtonBiomarker", {
     biomarker.vals <- reactiveValues(
      biomarker.list = suppressWarnings(find_biomarker(vals$MAE,
                    tax_level=input$taxl_biomarker,
                    input_select_target_biomarker=input$select_target_condition_biomarker,
                    nfolds = input$num.cv.nfolds,
                    nrepeats = input$num.biomarker.run,
                    seed = 99,
                    percent_top_biomarker = input$percent_top_biomarker,
                    model_name = input$select_model_biomarker))
     )


      output$biomarker_list <- renderTable({
        biomarker.vals$biomarker.list$biomarker
      })

      output$importance_plot <- renderPlot({
         biomarker.vals$biomarker.list$importance_plot
      })

      output$roc_plot <- renderPlot({
       biomarker.vals$biomarker.list$roc_plot
      })

  })
})
# 
# 
# 
#   ### check whether selected condition for biomarker has two levels or more
#   output$biomarker_condition_type <- reactive({
#     MAE = vals$MAE
#     target.var.index <- which(colnames(colData(MAE)) == input$select_target_condition_biomarker)
#     if (is_integer0(target.var.index)){
#       target.var.index <- 1
#     }
#     label.vec <- colData(MAE)[[target.var.index]]
#     label.level.num <- length(unique(label.vec))
#     if (label.level.num == 2){
#       return("binary")
#     } else{
#       return("multiple")
#     }
# 
#   })
#   outputOptions(output, "biomarker_condition_type", suspendWhenHidden = FALSE)
# 
#   # select 2 levels
#   output$biomarker_condition_options <- renderUI({
#     MAE = vals$MAE
#     variable.vec <- colData(MAE)[[
#       which(colnames(colData(MAE)) == input$select_target_condition_biomarker)]]
#     filter.option.vec <- sort(unique(variable.vec))
#     tagList(
#       selectInput("biomarker_condition_options_use", "Select 2 levels", choices = filter.option.vec, multiple = TRUE)
#     )
#   })
