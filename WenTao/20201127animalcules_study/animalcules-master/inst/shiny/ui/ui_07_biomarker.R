tabPanel("Biomarker",
  br(),
  sidebarLayout(
    sidebarPanel(
      br(),
      selectizeInput('taxl_biomarker', 'Taxonomy Level', choices = tax.name,
                     selected='genus'),
      selectInput("select_target_condition_biomarker", "Select Target Condition:",
                  covariates.two.levels),
        numericInput("num.cv.nfolds", "Number of CV nfolds", value = 3, max = 20, min = 3),
        numericInput("num.biomarker.run", "Number of CV repeats", value = 3, max = 100, min = 3),
        numericInput("percent_top_biomarker", "Top biomarker proportion", value = 0.2, max = 1, min = 0.01),
        selectInput("select_model_biomarker", "Select Model", c("logistic regression", "random forest")
      ),
      helpText("Note: we recommend to use this section only when sample size is larger than 100. 
                Smaller dataset might be biased to imbalanced class in split folds."),
      withBusyIndicatorUI(
        actionButton("goButtonBiomarker",
                     "Run",
                     class = "btn-primary")
      ),
      width=3
   ),
   mainPanel(
    tabsetPanel(
      tabPanel("Biomarker",
               br(),
               tableOutput("biomarker_list")
        ),
        tabPanel("Importance plot",
                 br(),
                 plotOutput("importance_plot",height=300)
        ),
        tabPanel("CV ROC plot",
                 br(),
                 plotOutput("roc_plot")
        )
      ), width=9
    )
  )
)
