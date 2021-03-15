tabPanel("Upload",
useShinyjs(),
tags$style(appCSS),
tags$div(
    class = "jumbotron",
    tags$div(
        class = "container",
        fluidRow(
            column(7, h1("animalcules"))
        ),
        p("Interactive Microbiome Analysis Toolkit"),
        uiOutput("tab")

    )
),
sidebarLayout(
    sidebarPanel(
        tags$span(style="color:#72bcd4", "Application Settings"),
        checkboxInput("global_adv", "Always Show Advanced"),
        checkboxInput("upload_adv", "Alternative Upload"),

        conditionalPanel(condition = "input.upload_adv == false",
            tags$span(style="color:#72bcd4", "Upload"),
            radioButtons("uploadChoice", "",
                        c("Count file" = "count",
                          "Example data" = "example",
                          "animalcules object" = "animalcules.file"
                        ))
            ),

       conditionalPanel(
            condition = "input.upload_adv == true",
            tags$span(style="color:#72bcd4", "Upload"),

            radioButtons("uploadChoiceAdv", "",
              c(
                "BIOM file" = "biom",
                "Count file with tax id" = 'countTi',
                "PathoScope file" = "pathofiles",
                "animalcules-id file" = "animalcules-id"
                
              ))
          ),
    conditionalPanel(
            condition = "input.upload_adv == false",
       conditionalPanel(condition = sprintf("input['%s'] == 'example'", "uploadChoice"),
                          selectInput("example_data", "Example dataset",
              c(
                "Simulated dataset" = "toy",
                "TB 16S profiling" = "tb",
                "Asthma metatranscriptomics" = "asthma"

                )),
                        withBusyIndicatorUI(
                          actionButton("upload_example",
                                       "Upload",
                                       class = "btn-primary")
                        )
       ),
       conditionalPanel(condition = sprintf("input['%s'] == 'animalcules.file'", "uploadChoice"),
                        fileInput("rdfile", ".rds file (required):",
                                  accept = c(
                                    ".rds"
                                  )
                        ),
                        radioButtons("rdtype", "Filetype",
                                     choices = c(
                                                 rds = "rds"
                                     ),
                                     selected = "rds"
                        ),
                        withBusyIndicatorUI(
                          actionButton("upload_animalcules",
                                       "Upload",
                                       class = "btn-primary")
                        )

       ),

       conditionalPanel(condition = sprintf("input['%s'] == 'count'", "uploadChoice"),
                        fileInput("countsfile", "Counts file (required):",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        fileInput("taxon.table", "Taxonomy table file (required):",
                                  accept = c(
                                    "text/csv",
                                    "text/comma-separated-values",
                                    "text/tab-separated-values",
                                    "text/plain",
                                    ".csv",
                                    ".tsv"
                                  )
                        ),
                        fileInput("annotfile.count", "Annotation file (required):",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        numericInput("metadata_sample_name_col_count", "Which column in metadata is sample name?",
                                     value = 1),
                        # Input: Checkbox if file has header ----
                        checkboxInput("header.count", "Header", TRUE),

                        # Input: Select separator ----
                        radioButtons("sep.count", "Separator",
                                     choices = c(Tab = "\t",
                                                 Comma = ",",
                                                 Semicolon = ";"
                                     ),
                                     selected = ","),
                        withBusyIndicatorUI(
                            actionButton("uploadDataCount",
                                         "Upload",
                                         class = "btn-primary")
                        )
       )
       ),
        conditionalPanel(
            condition = "input.upload_adv == true",
       conditionalPanel(condition = sprintf("input['%s'] == 'animalcules-id'", "uploadChoiceAdv"),
                        fileInput("rdfile_id", ".rds file (required):",
                                  accept = c(
                                    ".rds"
                                  )
                        ),
                        radioButtons("mae_data_type", "Choose count type",
                                     choices = c(
                                                 "EM count" = "em",
                                                 "Best hit" = 'hit'
                                     )
                        ),
                        withBusyIndicatorUI(
                          actionButton("upload_mae",
                                       "Upload",
                                       class = "btn-primary")
                        )

       ),
       conditionalPanel(condition = sprintf("input['%s'] == 'biom'", "uploadChoiceAdv"),
                        fileInput("biom_id", "biom file (required):",
                                  accept = c(
                                    ".biom"
                                  )
                        ),
                        helpText('Make sure the .biom file has sample metadata included.'),
                        withBusyIndicatorUI(
                          actionButton("upload_biom",
                                       "Upload",
                                       class = "btn-primary")
                        )

       ),
       conditionalPanel(condition = sprintf("input['%s'] == 'countTi'", "uploadChoiceAdv"),
                        fileInput("countsfileTi", "Counts file (required):",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        fileInput("annotfile.countTi", "Annotation file (required):",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        numericInput("metadata_sample_name_col_countTi", "Which column in metadata is sample name?",
                                     value = 1),
                        # Input: Checkbox if file has header ----
                        checkboxInput("header.countTi", "Header", TRUE),

                        # Input: Select separator ----
                        radioButtons("sep.countTi", "Separator",
                                     choices = c(Tab = "\t",
                                                 Comma = ",",
                                                 Semicolon = ";"
                                     ),
                                     selected = ","),
                        withBusyIndicatorUI(
                            actionButton("uploadDataCountTi",
                                         "Upload",
                                         class = "btn-primary")
                        )
       ),
       conditionalPanel(condition = sprintf("input['%s'] == 'pathofiles'", "uploadChoiceAdv"),
                        h5("Upload PathoScope generated .tsv files:"),
                        fileInput("countsfile.pathoscope", "PathoScope outputs (required):",
                                  multiple = TRUE,
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        fileInput("annotfile.ps", "Annotation file (required):",
                                  accept = c(
                                      "text/csv",
                                      "text/comma-separated-values",
                                      "text/tab-separated-values",
                                      "text/plain",
                                      ".csv",
                                      ".tsv"
                                  )
                        ),
                        textInput("report_suffix", "Report suffix", value = "-sam-report.tsv"),
                        numericInput("metadata_sample_name_col", "Which column in metadata is sample name?",
                                     value = 1),
                        # Input: Checkbox if file has header ----
                        checkboxInput("header.ps", "Header", TRUE),

                        # Input: Select separator ----
                        radioButtons("sep.ps", "Separator",
                                     choices = c(Tab = "\t",
                                                 Comma = ",",
                                                 Semicolon = ";"
                                     ),
                                     selected = "\t"),
                        withBusyIndicatorUI(
                            actionButton("uploadDataPs",
                                         "Upload",
                                         class = "btn-primary")
                        ),
                        helpText("This might take 10-20 seconds to upload.")

       )
        )
   ),
   mainPanel(
             conditionalPanel(
            condition = "input.upload_adv == true",
       conditionalPanel(condition = "input.uploadChoiceAdv === 'pathofiles'",
                        h4("Note: please click \"open in browser\" for enabling functions like multiple files upload."),
                        helpText("Counts Table: column names must be sample name"),
                        DT::dataTableOutput("contents.count"),
                        helpText("Annotation table"),
                        DT::dataTableOutput("contents.meta")
       ),
       conditionalPanel(condition = "input.uploadChoiceAdv === 'countTi'",

                        #tags$img(src='count_table_example.png', height = 180, width = 800),
                        helpText("Counts Table: "),
                        helpText("1. Column names must be sample name"),
                        helpText("2. The first column must be Tax id"),

                        DT::dataTableOutput("contents.count.2Ti"),

                        helpText("Annotation table: "),
                        helpText("1. Row names must be sample name"),
                        helpText("2. The first row must sample attribute labels"),

                        DT::dataTableOutput("contents.meta.2Ti")
        ),
       
       conditionalPanel(condition = "input.uploadChoiceAdv === 'biom'",
                        h5('Note: Please check http://biom-format.org/documentation/adding_metadata.html 
                                 to add sample metadate into .biom file if sample metadata is missing.'),
                        helpText("Counts Table"),
                        DT::dataTableOutput("biom.count"),
                        helpText("Annotation table"),
                        DT::dataTableOutput("biom.meta"),
                        helpText("Taxonomy table"),
                        DT::dataTableOutput("biom.tax")
       )),
       
             conditionalPanel(
            condition = "input.upload_adv == false",
       conditionalPanel(condition = "input.uploadChoice === 'example'",
                        
              conditionalPanel(
                  condition = "input.example_data == 'tb'",
                  h4("TB 16S profiling"),
                  h5("TB 16S profiling is a real dataset containing 30 samples and 417 microbes. "),
                 h6("Reference: Botero LE, Delgado-Serrano L, Cepeda ML, Bustos JR, Anzola JM, 
                 Del Portillo P, Robledo J, Zambrano MM. Respiratory tract clinical sample selection for 
                 microbiota analysis in patients with pulmonary tuberculosis. Microbiome. 2014 Aug 25. 
                 doi: 10.1186/2049-2618-2-29.")
              ),
              conditionalPanel(
                  condition = "input.example_data == 'asthma'",
                  h4("Asthma metatranscriptomics"),
                  h5("Nasal swabs metatranscriptomic data from 8 asthma children and 6 control children."),
                 h6("Castro-Nallar E, Shen Y, Freishtat RJ, Pérez-Losada M, Manimaran S, Liu G, 
                Johnson WE, Crandall KA. Integrating microbial and host transcriptomics to 
                characterize asthma-associated microbial communities. BMC Med Genomics. 
                2015 Aug 16;8:50. doi: 10.1186/s12920-015-0121-1. PubMed PMID: 26277095; 
                PubMed Central PMCID: PMC4537781.")
              ),
              conditionalPanel(
                  condition = "input.example_data == 'toy'",
                  h4("Simulated dataset"),
                  h5("Simulated dataset is a small synthetic microbiome dataset containing 50 samples and 100 microbes. 
                     It's loaded already, so you could simply continue and play with animalcules from this point!")
                  
              )
  
       ),
       conditionalPanel(condition = "input.uploadChoice === 'count'",

                        tags$img(src='count_table_example.png', height = 180, width = 800),
                        helpText("Counts Table: "),
                        helpText("1. Column names must be sample name"),
                        helpText("2. The first column must be microbe name"),

                        DT::dataTableOutput("contents.count.2"),
                        helpText("Taxonomy Table: "),
                        helpText("1. Column names must be taxonomy levels, like family, genus, species..."),
                        helpText("2. The first column must be microbe name"),

                        DT::dataTableOutput("contents.taxonomy"),
                        helpText("Annotation table: "),
                        helpText("1. Row names must be sample name"),
                        helpText("2. The first row must sample attribute labels"),

                        DT::dataTableOutput("contents.meta.2")
        )
      ))
   )

)
