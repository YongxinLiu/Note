tabPanel("Summary and Filter",
  tabsetPanel(
    tabPanel("Filter",
      br(),
      sidebarLayout(
        sidebarPanel(
          # View Style
          selectInput("filter_type", "Filter By", c("Metadata", "Microbes")),

          # Metadata
          conditionalPanel(condition = "input.filter_type == 'Metadata'",
            selectizeInput('filter_type_metadata', 'Select a Condition', choices=covariates, multiple=FALSE),
            uiOutput("filter_metadata_params"),
            withBusyIndicatorUI(
              actionButton("filter_metadata_btn", "Filter", class = "btn-primary")
            )
          ),

          # Microbes
          conditionalPanel(condition = "input.filter_type == 'Microbes'",
            selectInput("filter_type_microbes", "Select Filter Condition", c("Average Read Number",
                                                                             "Average Relative Abundance",
                                                                             "Average Prevalence")
            )
          ),
          conditionalPanel(condition = "input.filter_type == 'Microbes' & input.filter_type_microbes == 'Average Read Number'",
            numericInput("filter_microbes_read_inp", "Set Minimum", 0, min = 0, max = 10000),
            withBusyIndicatorUI(
              actionButton("filter_microbes_read_btn", "Filter", class = "btn-primary")
            )
          ),
          conditionalPanel(condition = "input.filter_type == 'Microbes' & input.filter_type_microbes == 'Average Relative Abundance'",
            numericInput("filter_microbes_rela_min_inp", "Set Minimum", 0, min = 0, max = 1),
            numericInput("filter_microbes_rela_max_inp", "Set Maximum", 1, min = 0, max = 1),
            withBusyIndicatorUI(
              actionButton("filter_microbes_rela_btn", "Filter", class = "btn-primary")
            )
          ),
          conditionalPanel(condition = "input.filter_type == 'Microbes' & input.filter_type_microbes == 'Average Prevalence'",
            numericInput("filter_microbes_prev_min_inp", "Set Minimum", 0, min = 0, max = 1),
            numericInput("filter_microbes_prev_max_inp", "Set Maximum", 1, min = 0, max = 1),
            withBusyIndicatorUI(
              actionButton("filter_microbes_prev_btn", "Filter", class = "btn-primary")
            )
          ),

          br(),

          # Reset
          withBusyIndicatorUI(
            actionButton("filter_reset_btn", "Reset")
          ),
          width=5,
          br(),
          downloadButton('download_rds', 'Download Animalcules File'),
          br(),
          downloadButton('download_biom', 'Download Biom File'),


          checkboxInput("filter_adv", "Advanced Options"),
          conditionalPanel(
            condition = "input.filter_adv == true | input.global_adv == true",
          # Discard samples
          selectizeInput("filter_sample_dis", "Discard Samples", choices=sam.name, multiple=TRUE),
          # Discard organisms
          selectizeInput("filter_organism_dis", "Discard Organisms", choices=org.name, multiple=TRUE),
          withBusyIndicatorUI(
            actionButton("filter_discard_btn", "Discard")
          )
          ),


          br()

        ),
        mainPanel(
          fluidRow(
            column(5,
              tableOutput("filter_summary_table")
            ),
            column(7,
              plotlyOutput("filter_summary_top_plot", height="350px"),
              plotlyOutput("filter_summary_bottom_plot", height="350px")
            )
          ),
          width=7
        )
      )
    ),
    tabPanel("Categorize",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput('filter_bin_cov', 'Continuous Variables', choices=num_covariates, multiple=FALSE),
          uiOutput("filter_nbins"),

          textInput("filter_new_covariate", "New Variable Label", value = "new_var"),

          checkboxInput("filter_bin_adv", "Advanced Options"),

          conditionalPanel(
            condition = "input.filter_bin_adv == true | input.global_adv == true",
            textInput('filter_bin_breaks', 'Custom Breaks (Comma Delimited)'),
            verbatimTextOutput("filter_bin_to1"),
            textInput('filter_bin_labels', 'Custom Labels (Comma Delimited)'),
            verbatimTextOutput("filter_bin_to2")
          ),

          actionButton("filter_create_bins", "Create Bins", class = "btn-primary"),
          width=5
        ),
        mainPanel(
          plotlyOutput("filter_unbin_plot", height="200px"),
          br(),
          plotlyOutput("filter_bin_plot"),
          width=7
        )
      )
    ),
    tabPanel("Assay Dashboard",
      br(),
      sidebarLayout(
        sidebarPanel(
          selectInput("select_assay", "Select Assay",
              c("Count table" = "count",
                "Relative Abundance table" = "ra",
                "LogCPM table" = "logcpm",
                "Taxonomy table" = "tax",
                "Annotation table" = "annot"
              )),
        conditionalPanel(condition = sprintf("input['%s'] == 'count'", "select_assay"),
          selectizeInput("assay_count_taxlev", "Tax Level", choices=tax.name, selected=tax.default),
          actionButton("view_assay_count",
             "View",
             class = "btn-primary"),
          downloadButton('download_assay_count', 'Download')
        ),
        
        conditionalPanel(condition = sprintf("input['%s'] == 'ra'", "select_assay"),
          selectizeInput("assay_ra_taxlev", "Tax Level", choices=tax.name, selected=tax.default),
          actionButton("view_assay_ra",
             "View",
             class = "btn-primary"),
          downloadButton('download_assay_ra', 'Download')
        ),
        conditionalPanel(condition = sprintf("input['%s'] == 'logcpm'", "select_assay"),
          selectizeInput("assay_logcpm_taxlev", "Tax Level", choices=tax.name, selected=tax.default),
          actionButton("view_assay_logcpm",
             "View",
             class = "btn-primary"),
          downloadButton('download_assay_logcpm', 'Download')
        ),
        conditionalPanel(condition = sprintf("input['%s'] == 'tax'", "select_assay"),
          actionButton("view_assay_tax",
             "View",
             class = "btn-primary"),
          downloadButton('download_assay_tax', 'Download')
        ),
        conditionalPanel(condition = sprintf("input['%s'] == 'annot'", "select_assay"),
          actionButton("view_assay_annot",
             "View",
             class = "btn-primary"),
          downloadButton('download_assay_annot', 'Download')
        )
        
        ),
        mainPanel(
          conditionalPanel(condition = sprintf("input['%s'] == 'count'", "select_assay"),
              DT::dataTableOutput("assay_table_count", width='95%')
        
          ),
          conditionalPanel(condition = sprintf("input['%s'] == 'ra'", "select_assay"),
              DT::dataTableOutput("assay_table_ra", width='95%')
        
          ),
          conditionalPanel(condition = sprintf("input['%s'] == 'logcpm'", "select_assay"),
              DT::dataTableOutput("assay_table_logcpm", width='95%')
        
          ),
          conditionalPanel(condition = sprintf("input['%s'] == 'tax'", "select_assay"),
              DT::dataTableOutput("assay_table_tax", width='95%')
        
          ),
          conditionalPanel(condition = sprintf("input['%s'] == 'annot'", "select_assay"),
              DT::dataTableOutput("assay_table_annot", width='95%')
        
          )
          
          
        )
      )
    )
  )
)
