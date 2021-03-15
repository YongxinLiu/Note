tabPanel("Dimension Reduction",
  tabsetPanel(
    tabPanel("PCA",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          selectizeInput('dimred_pca_taxlev', 'Taxonomy Level', choices = tax.name, selected=tax.default),

          selectInput("dimred_pca_color", "Color points by:", covariates),

          checkboxInput("dimred_pca_adv", "Advanced Options (3D)"),

          conditionalPanel(
            condition = "input.dimred_pca_adv == true | input.global_adv == true",
            numericInput('dimred_pca_x', 'Principal Component (x-axis)', 1, min=1, max=50),
            numericInput('dimred_pca_y', 'Principal Component (y-axis)', 2, min=1, max=50),
            numericInput('dimred_pca_z', 'Principal Component (z-axis)', NA, min=1, max=50),
            selectInput("dimred_pca_shape", "Shape points by:", c("None", covariates.colorbar)),
            selectInput("dimred_pca_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                     "Counts"             = "counts",
                                                                     "log(CPM)"           = "logcpm"),
                                                                     selected             = "logcpm")
          ),

          # Do plot button
          actionButton("dimred_pca_plot_btn", "Plot", class = "btn-primary"),
          actionButton("dimred_pca_table_btn", "Table"),
          width=3
        ),
        mainPanel(
          fluidRow(
            column(7,
              plotlyOutput("dimred_pca_plot", height="500px")
            ),
            column(5,
              dataTableOutput("dimred_pca_table")
            )
          ),
          width=9
        )
      )
    ),
    tabPanel("PCoA",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          selectizeInput('dimred_pcoa_taxlev', 'Taxonomy Level', choices = tax.name, selected=tax.default),

          selectInput("dimred_pcoa_color", "Color points by:", covariates),

          checkboxInput("dimred_pcoa_adv", "Advanced Options (3D)"),

          conditionalPanel(
            condition = "input.dimred_pcoa_adv == true | input.global_adv == true",
            numericInput('dimred_pcoa_x', 'Principal Coordinate (x-axis)', 1, min=1, max=50),
            numericInput('dimred_pcoa_y', 'Principal Coordinate (y-axis)', 2, min=1, max=50),
            numericInput('dimred_pcoa_z', 'Principal Coordinate (z-axis)', NA, min=1, max=50),
            selectInput("dimred_pcoa_shape", "Shape points by:", c("None", covariates.colorbar)),
            selectInput("dimred_pcoa_method", "Select distance metric",
                        c("bray", "jaccard"), selected = "bray")
          ),

          # Do plot button
          actionButton("dimred_pcoa_plot_btn", "Plot", class = "btn-primary"),
          actionButton("dimred_pcoa_table_btn", "Table"),
          width=3
        ),
        mainPanel(
          fluidRow(
            column(7,
              plotlyOutput("dimred_pcoa_plot", height="500px")
            ),
            column(5,
              dataTableOutput("dimred_pcoa_table")
            )
          ),
          width=9
        )
      )
    ),
    tabPanel("UMAP",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          selectizeInput('dimred_umap_taxlev', 'Taxonomy Level', choices = tax.name, selected=tax.default),

          selectInput("dimred_umap_color", "Color points by:", covariates),

          checkboxInput("dimred_umap_adv", "Advanced Options (3D)"),

          conditionalPanel(
            condition = "input.dimred_umap_adv == true | input.global_adv == true",
            numericInput('dimred_umap_x', 'Component (x-axis)', 1, min=1, max=50),
            numericInput('dimred_umap_y', 'Component (y-axis)', 2, min=1, max=50),
            numericInput('dimred_umap_z', 'Component (z-axis)', NA, min=1, max=50),
            numericInput('dimred_umap_n_neighbors', 'Nearest Neighbors', 15),
            selectizeInput('dimred_umap_metric', 'Distance Metric', c("euclidean", "manhattan"), selected="euclidean"),
            numericInput('dimred_umap_n_epochs', 'Iterations', 200),
            selectizeInput('dimred_umap_init', 'Initial Embedding', c("spectral", "random"), selected="spectral"),
            numericInput('dimred_umap_min_dist', 'Min Distance', 0.1),
            selectInput("dimred_umap_shape", "Shape points by:", c("None", covariates.colorbar)),
            selectInput("dimred_umap_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                      "Counts"             = "counts",
                                                                      "log(CPM)"           = "logcpm"),
                                                                      selected             = "logcpm")
          ),

          # Do plot button
          actionButton("dimred_umap_plot_btn", "Plot", class = "btn-primary"),
          width=3
        ),
        mainPanel(
          plotlyOutput("dimred_umap_plot", width="700px", height="700px"),
          width=9
        )
      )
    ),
    tabPanel("t-SNE",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          selectizeInput('dimred_tsne_taxlev', 'Taxonomy Level', choices = tax.name, selected=tax.default),

          selectInput("dimred_tsne_color", "Color points by:", covariates),

          checkboxInput("dimred_tsne_adv", "Advanced Options (3D)"),
          conditionalPanel(
            condition = "input.dimred_tsne_adv == true | input.global_adv == true",
            checkboxInput("dimred_tsne_cached", "Use Cached Data", TRUE),
            helpText("Note: Uncheck the \"Use Cached Data\" first to run 3D tSNE"),
            selectInput("dimred_tsne_k", "Select Final Dimensions", c("2D" = "2D",
                                                                      "3D" = "3D"),
                                                                      selected = "2D"),
            selectInput("dimred_tsne_shape", "Shape points by:", c("None", covariates.colorbar)),
            helpText("Note: Uncheck the \"Use Cached Data\" first to change perplexity"),
            numericInput('dimred_tsne_perplexity', 'Perplexity', 10, min=1, max=1000, step=0.1),
            numericInput('dimred_tsne_initial_dims', 'Inital Dimensions', 30, min=10, max=1000),
            selectInput("dimred_tsne_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                     "Counts"              = "counts",
                                                                     "log(CPM)"            = "logcpm"),
                                                                     selected              = "logcpm")
          ),

          # Do plot button
          withBusyIndicatorUI(
          actionButton("dimred_tsne_plot_btn", "Plot",class = "btn-primary")
          ),
          width=3
        ),
        mainPanel(
          helpText("Note: it will take longer time for the first run. 
         Subsequent plots will use cached assay unless the
         \"Use Cached Data\" is disabled from the advanced option."),
          fluidRow(
            plotlyOutput("dimred_tsne_plot", width="700px", height="700px")
          ),
          width=9
        )
      )
    )
  )
)
