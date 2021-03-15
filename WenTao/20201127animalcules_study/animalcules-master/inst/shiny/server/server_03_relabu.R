#
# Relative Abundance Barplot
#
# Plot when button is pressed
do_relabu_bar <- eventReactive(input$relabu_bar_plot_btn, {
    withBusyIndicatorServer("relabu_bar_plot_btn", {
    p <- relabu_barplot(MAE = vals$MAE,
                        tax_level = input$relabu_bar_taxlev,
                        order_organisms = input$relabu_bar_org_order,
                        sort_by = input$relabu_bar_sort,
                        group_samples = input$relabu_bar_group_samples,
                        group_conditions = input$relabu_bar_group_conditions,
                        sample_conditions = input$relabu_bar_sample_conditions,
                        isolate_samples = input$relabu_bar_sample_iso,
                        discard_samples = input$relabu_bar_sample_dis,
                        show_legend = input$relabu_bar_legend)
    return(p)
    })
})

# Reaction to button pressing
  output$relabu_bar_plot <- renderPlotly({
  p <- do_relabu_bar()
    return(p)
})

# Used to dynamically adjust plot height/width
output$relabu_bar_dynamic_plot <- renderUI({
    height = paste(input$relabu_bar_height, "px", sep="")
    width = paste(input$relabu_bar_width, "px", sep="")
    plotlyOutput("relabu_bar_plot", width=width, height=height)
})

# Return unique organisms for a given tax level
output$relabu_bar_org_order <- renderUI({
    organisms <- unique(as.data.frame(rowData(experiments(vals$MAE)[[1]]))[,input$relabu_bar_taxlev])
    selectizeInput('relabu_bar_org_order', label='Order Organisms', choices=organisms, multiple=TRUE)
})

#
# Relative Abundance Heatmap
#
# Plot when button is pressed
do_relabu_heatmap <- eventReactive(input$relabu_heatmap_plot_btn, {
    withBusyIndicatorServer("relabu_heatmap_plot_btn", {
    p <- relabu_heatmap(MAE = vals$MAE,
                        tax_level = input$relabu_heatmap_taxlev,
                        sort_by = input$relabu_heatmap_sort,
                        sample_conditions = input$relabu_heatmap_conditions,
                        isolate_organisms = input$relabu_heatmap_org_iso,
                        isolate_samples = input$relabu_heatmap_sample_iso,
                        discard_samples = input$relabu_heatmap_sample_dis,
                        log_cpm = input$relabu_heatmap_logcpm)
    return(p)
    })
})

# Reaction to button pressing
output$relabu_heatmap_plot <- renderPlotly({
    p <- do_relabu_heatmap()
    return(p)
})

# Used to dynamically adjust plot height/width
output$relabu_heatmap_dynamic_plot <- renderUI({
    height = paste(input$relabu_heatmap_height, "px", sep="")
    width = paste(input$relabu_heatmap_width, "px", sep="")
    plotlyOutput("relabu_heatmap_plot", width=width, height=height)
})

# Return unique organisms for a given tax level
output$relabu_heatmap_org_iso <- renderUI({
    organisms <- unique(as.data.frame(rowData(experiments(vals$MAE)[[1]]))[,input$relabu_heatmap_taxlev])
    selectizeInput('relabu_heatmap_org_iso', label='Isolate Organisms', choices=organisms, multiple=TRUE)
})

#
# Relative Abundance Boxplot
#
# Plot when button is pressed
do_relabu_box <- eventReactive(input$relabu_box_plot_btn, {
    withBusyIndicatorServer("relabu_box_plot_btn", {
    tavlevs <- as.list(input$relabu_box_taxlev)
    plots <- lapply(tavlevs, function(x) {
        id <- paste("relabu_box_organisms", x, sep="_")
        organisms <- input[[id]]

        # One organism one plot
        if (length(organisms) == 1) {
            relabu_boxplot(MAE = vals$MAE,
                           tax_level = x,
                           condition = input$relabu_box_condition,
                           organisms = organisms,
                           datatype = input$relabu_box_datatype) %>%
            # Dynamic sizing
            plotly::layout(height=input$relabu_box_height, width=input$relabu_box_width)
        }
        # Multiple organisms multiple plots
        else if (length(organisms) > 1) {
            # Merged plots
            if (!input$relabu_box_separate) {
                relabu_boxplot(MAE = vals$MAE,
                               tax_level = x,
                               condition = input$relabu_box_condition,
                               organisms = organisms,
                               datatype = input$relabu_box_datatype) %>%
                # Dynamic sizing
                plotly::layout(height=input$relabu_box_height, width=input$relabu_box_width)
            } 
            # Separate plots
            else {
                subplots <- lapply(as.list(organisms), function(y) {
                    relabu_boxplot(MAE = vals$MAE,
                                   tax_level = x,
                                   condition = input$relabu_box_condition,
                                   organisms = y,
                                   datatype = input$relabu_box_datatype)
                })
                # Keep only one legend
                for (i in 2:length(subplots)) {
                    subplots[[i]] %<>% style(showlegend = FALSE)

                }
                # Subplot and Dynamic sizing
                plotly::subplot(subplots) %>% layout(height=input$relabu_box_height, width=input$relabu_box_width)
            }
        }
        # No organisms selected
        else {
            plotly_empty()
        }
    })
    return(plots)
    })
})

# Reaction to button pressing
output$relabu_box_plots <- renderUI({
    p <- do_relabu_box()
    return(p)
})

# Return a dynamic number of organism choices for each tax level selected
output$relabu_box_organisms <- renderUI({
    tavlevs <- as.list(input$relabu_box_taxlev)
    inputs <- lapply(tavlevs, function(x) {
        id <- paste("relabu_box_organisms", x, sep="_")
        organisms <- unique(as.data.frame(rowData(experiments(vals$MAE)[[1]]))[,x])
        selectizeInput(id, label=x, choices=organisms, selected=organisms[1], multiple=TRUE)
    })
    return(inputs)
})
