tabPanel("Abundance",
  tabsetPanel(
    tabPanel("Barplots",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          # Sort the samples by a condition
          conditionalPanel(
            condition = "input.relabu_bar_group_samples == false",
            selectizeInput('relabu_bar_sample_conditions', 'Color Samples by Condition', choices=covariates, multiple=TRUE)
          ),
          conditionalPanel(
            condition = "input.relabu_bar_group_samples == true",
            selectizeInput('relabu_bar_group_conditions', 'Color Samples by Condition', choices=c("ALL", covariates))
          ),

          # Sample aggregation
          checkboxInput("relabu_bar_group_samples", "Group Samples by Condition"),

          # Select taxon level
          selectInput("relabu_bar_taxlev", "Tax Level", choices=tax.name, selected=tax.default),

          # Sort the bars
          radioButtons("relabu_bar_sort", "Sort By", c("No Sorting" = "nosort",
                                                       "Conditions" = "conditions",
                                                       "Organisms" = "organisms",
                                                       "Alphabetically" = "alphabetically"),
                                                        selected = "nosort"),

          # Advanced options
          checkboxInput("relabu_bar_adv", "Advanced Options"),

          # Isolate samples
          conditionalPanel(
            condition = "input.relabu_bar_adv == true | input.global_adv == true",
            selectizeInput("relabu_bar_sample_iso", "Isolate Samples", choices=sam.name, multiple=TRUE),
            selectizeInput("relabu_bar_sample_dis", "Discard Samples", choices=sam.name, multiple=TRUE),
            uiOutput("relabu_bar_org_order"),
            checkboxInput("relabu_bar_legend", "Show Legend", value=TRUE),
            sliderInput("relabu_bar_height", "Plot Height", 400, 1200, value=600, step=50, post="px"),
            sliderInput("relabu_bar_width", "Plot Width", 400, 1200, value=800, step=50, post="px")
          ),

          # Do plot button
          withBusyIndicatorUI(
          actionButton("relabu_bar_plot_btn", "Plot", class = "btn-primary")
          ),
          width=3
        ),
        mainPanel(
          uiOutput("relabu_bar_dynamic_plot"),
          width=9
        )
      )
    ),
    tabPanel("Heatmap",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          # Sort the samples by a condition
          selectizeInput("relabu_heatmap_conditions", "Color Samples by Condition", choices=covariates, multiple=TRUE),

          # Select taxon level
          selectInput("relabu_heatmap_taxlev", "Tax Level", choices=tax.name, selected=tax.default),


          # Column sort
          radioButtons("relabu_heatmap_sort", "Sort By", c("No Sorting" = "nosort",
                                                           "Conditions" = "conditions",
                                                           "Organisms" = "organisms",
                                                           "Alphabetically" = "alphabetically"),
                                                           selected = "nosort"),

          # Advanced options
          checkboxInput("relabu_heatmap_adv", "Advanced Options"),

          # Dynamically generate based on tax level
          conditionalPanel(
            condition = "input.relabu_heatmap_adv == true | input.global_adv == true",
            uiOutput("relabu_heatmap_org_iso"),
            selectizeInput("relabu_heatmap_sample_iso", "Isolate Samples", choices=sam.name, multiple=TRUE),
            selectizeInput("relabu_heatmap_sample_dis", "Discard Samples", choices=sam.name, multiple=TRUE),
            checkboxInput("relabu_heatmap_logcpm", "log(CPM)", value=TRUE),
            sliderInput("relabu_heatmap_height", "Plot Height", 400, 1200, value=600, step=50, post="px"),
            sliderInput("relabu_heatmap_width", "Plot Width", 600, 1000, value=800, step=50, post="px")
          ),

          # Do plot button
          withBusyIndicatorUI(
          actionButton("relabu_heatmap_plot_btn", "Plot", class = "btn-primary")
          ),
          width=3
        ),
        mainPanel(
          uiOutput("relabu_heatmap_dynamic_plot"),
          width=9
        )
      )
    ),
    tabPanel("Boxplots",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          selectInput("relabu_box_taxlev", "Tax Levels", choices=tax.name, selected=tax.default, multiple=FALSE),

          # Dynamic choose from organisms based on tax level
          uiOutput("relabu_box_organisms"),

          # Separate plots
          checkboxInput("relabu_box_separate", "Separate Plots"),

          # Select condition
          selectInput("relabu_box_condition", "Select condition", covariates.colorbar),

          # Select datatype
          radioButtons("relabu_box_datatype", "Select data format", c("Relative Abundance" = "relative abundance",
                                                                      "Counts"             = "counts",
                                                                      "log(CPM)"           = "logcpm"),
                                                                      selected             = "logcpm"),

          # Advanced options
          checkboxInput("relabu_box_adv", "Advanced Options"),

          # Adjust height of plot
          conditionalPanel(
            condition = "input.relabu_box_adv == true | input.global_adv == true",
            sliderInput("relabu_box_height", "Plot Height", 200, 1000, value=400, step=50, post="px"),
            sliderInput("relabu_box_width", "Plot Width", 600, 1400, value=1000, step=50, post="px")
          ),

          # Do plot button
          withBusyIndicatorUI(
            actionButton("relabu_box_plot_btn", "Plot", class = "btn-primary")
          ),
          width=3
        ),
        mainPanel(
          uiOutput("relabu_box_plots"),
          width=9
        )
      )
    )
  )
)
