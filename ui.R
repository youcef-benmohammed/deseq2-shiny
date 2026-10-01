navbarPage(
  title = div(icon("dna"), "DESeq2 Shiny"),
  id = "main_nav",
  theme = shinytheme("flatly"),
  header = tagList(
    useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", href = "style.css"))
  ),
  tabPanel(div(icon("upload"), "Data & model"), mod_data_ui("data")),
  tabPanel(div(icon("chart-bar"), "QC"), mod_qc_ui("qc")),
  tabPanel(div(icon("table"), "Normalised counts"), mod_norm_ui("norm")),
  tabPanel(div(icon("filter"), "Differential expression"), mod_de_ui("de"))
)
