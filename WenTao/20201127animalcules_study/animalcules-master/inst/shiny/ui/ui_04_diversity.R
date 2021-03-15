tabPanel("Diversity",
  tabsetPanel(
    tabPanel("Alpha Diversity",
      br(),
      sidebarLayout(
        sidebarPanel(
          br(),
          selectizeInput('taxl.alpha', 'Taxonomy Level', choices = tax.name, selected=tax.default),
          selectInput("select_alpha_div_condition", "Compare between:", covariates),
          selectInput("select_alpha_div_method", "Choose method:", alpha.methods),
          actionButton("alpha_boxplot", "Run", class = "btn-primary")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Boxplot",
              plotlyOutput("AlphaDiversity"),
              br(),
              DT::dataTableOutput("alpha.stat.test")
            )
            # ,tabPanel("Alpha Diversity Table",
            #   br(),
            #   DT::dataTableOutput("table.alpha"),
            #   downloadButton('download_table_alpha', 'Download')
            # )
          )
        )
      )
    ),
    tabPanel("Beta Diversity",
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput('taxl.beta', 'Taxonomy Level', choices = tax.name, selected=tax.default),
          selectizeInput('bdhm_select_conditions', 'Color Samples by Condition', choices=covariates.colorbar, multiple=TRUE),
          radioButtons("bdhm_sort_by", "Sort By", c("No Sorting" = "nosort", "Conditions" = "conditions"), selected="nosort"),
          selectInput("beta_method", "Choose distance metric:", choices = c("bray", "jaccard", 'unweighted unifrac', 'weighted unifrac'), selected="bray"),
          actionButton("beta_heatmap", "Plot Heatmap", class = "btn-primary")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Heatmap",
              br(),
              plotlyOutput("BetaDiversityHeatmap", width="800px", height="600px")
            ),
            tabPanel("Boxplot",
              fluidRow(
                column(4,
                  selectInput("select_beta_stat_method","Select Test", c("PERMANOVA", "Kruskal-Wallis", "Wilcoxon rank sum test")),
                  helpText("Only variables with 2 levels are supported" ),
                  selectInput("select_beta_condition", "Select condition", covariates.two.levels),
                  numericInput("num_permutation_permanova", "Number of permutations", value = 999, max = 2000)
                ),
                column(8,
                  br(),
                  DT::dataTableOutput("beta.stat.test"),
                  actionButton("beta_boxplot", "Run", class = "btn-primary")
                )
              ),
              plotlyOutput("BetaDiversityBoxplot")
            )
            # ,tabPanel("Beta Diversity Table",
            #   br(),
            #   DT::dataTableOutput("table.beta"),
            #   downloadButton('download_table_beta', 'Download')
            # )
          )
        )
      )
    )
  )
)
