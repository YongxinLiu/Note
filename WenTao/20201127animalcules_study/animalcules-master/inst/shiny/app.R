library(shiny)
library(shinyjs)
library(MultiAssayExperiment)
library(ggplot2)
library(plotly)
library(vegan)
library(dplyr)
library(magrittr)
library(biomformat)

source(file.path("utils", "helpers.R"),  local = TRUE)
# source(file.path("utils", "server_util.R"),  local = TRUE)
source(file.path("utils", "ui_util.R"),  local = TRUE)

ui <- navbarPage(
  title = paste("animalcules v", packageVersion("animalcules"), sep = ""),
  id="Animalcules",
  fluid=TRUE,
  theme = "bootstrap.min.css",
  source(file.path("ui", "ui_01_upload.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_02_filter.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_03_relabu.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_04_diversity.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_05_dimred.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_06_differential.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_07_biomarker.R"),  local = TRUE)$value
)

server <- function(input, output, session) {
  source(file.path("server", "server_01_upload.R"),  local = TRUE)$value
  source(file.path("server", "server_02_filter.R"),  local = TRUE)$value
  source(file.path("server", "server_03_relabu.R"),  local = TRUE)$value
  source(file.path("server", "server_04_diversity.R"),  local = TRUE)$value
  source(file.path("server", "server_05_dimred.R"),  local = TRUE)$value
  source(file.path("server", "server_06_differential.R"),  local = TRUE)$value
  source(file.path("server", "server_07_biomarker.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)
