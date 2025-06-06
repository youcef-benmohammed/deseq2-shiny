library(shiny)
library(DT)
library(shinyjs)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
library(shinythemes)

shinyUI(
  navbarPage(
    title = div(icon("dna"), "DESeq2 Analysis"),
    theme = shinytheme("flatly"),
    useShinyjs(),
    
    tabPanel(
      title = div(icon("upload"), "Upload Counts"),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "counts_file",
            tagList(icon("table"), "Count Matrix (TSV/TXT)"),
            accept = c(".tsv", ".txt")
          ),
          fileInput(
            "design_matrix",
            tagList(icon("project-diagram"), "Design Matrix (TSV/TXT)"),
            accept = c(".tsv", ".txt")
          ),
          hr(),
          actionButton(
            inputId = "run",
            label = tagList(icon("play"), "Run DESeq2"),
            class = "btn-primary",
            style = "width: 100%; margin-top: 10px;"
          ),
          hr(),
          div(id = "log", style = "font-size: 0.9em; color: #555; max-height: 150px; overflow-y: auto;"),
          width = 4
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Counts Table", DT::dataTableOutput("counts_table")),
            tabPanel("Design Table", DT::dataTableOutput("design_matrix_table"))
          ),
          width = 8
        )
      )
    ),
    
    tabPanel(
      title = div(icon("chart-bar"), "QC Plots"),
      sidebarLayout(
        sidebarPanel(
          uiOutput("group_by_select"),
          uiOutput("filter_col_select"),
          uiOutput("filter_values_ui"),
          uiOutput("pc_x_select"),
          uiOutput("pc_y_select"),
          width = 3
        ),
        mainPanel(
          plotlyOutput("pcaPlot", height = "400px"),
          hr(),
          plotlyOutput("elbowPlot", height = "300px"),
          hr(),
          plotlyOutput("sampleClustering", height = "600px"),
          width = 9
        )
      )
    ),
    
    tabPanel(
      title = div(icon("table"), "Normalized Counts"),
      sidebarLayout(
        sidebarPanel(
          downloadButton("normcounts", "Download Normalized Counts", class = "btn-success"),
          width = 3
        ),
        mainPanel(
          DT::dataTableOutput("normalized_counts"),
          width = 9
        )
      )
    ),
    
    tabPanel(
      title = div(icon("filter"), "Differential Expression"),
      sidebarLayout(
        sidebarPanel(
          uiOutput("condition_a"),
          uiOutput("condition_b"),
          hr(),
          sliderInput(
            "padj_thresh",
            "Adjusted p-value threshold",
            min = 0,
            max = 0.1,
            value = 0.05,
            step = 0.001
          ),
          sliderInput(
            "logfc_thresh",
            "Log2 Fold Change threshold",
            min = 0,
            max = 5,
            value = 1,
            step = 0.5
          ),
          hr(),
          actionButton("run_deg", tagList(icon("magic"), "Run DEG Analysis")),
          width = 3
        ),
        mainPanel(
          verbatimTextOutput("degSummary"),
          hr(),
          plotlyOutput("MAPlot", height = "400px"),
          hr(),
          plotlyOutput("pvalues", height = "300px"),
          hr(),
          plotlyOutput("volcanoPlot", height = "400px"),
          hr(),
          downloadButton("downloadDEG", "Download DEG Table"),
          br(), br(),
          DT::dataTableOutput("degTable"),
          width = 9
        )
      )
    )
  )
)
