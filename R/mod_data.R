# Module: data upload, validation and DESeq2 model fitting.

mod_data_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4(icon("upload"), "1. Input data"),
      fileInput(ns("counts_file"), "Count matrix (TSV/TXT)", accept = c(".tsv", ".txt")),
      fileInput(ns("design_file"), "Design matrix (TSV/TXT)", accept = c(".tsv", ".txt")),
      actionButton(ns("load_example"), tagList(icon("flask"), "Load example dataset (pasilla)"),
                   class = "btn-default", style = "width: 100%;"),
      helpText("Example: Drosophila pasilla knock-down (Brooks et al., 2011), 7 samples."),
      hr(),
      h4(icon("sliders-h"), "2. Model"),
      uiOutput(ns("model_ui")),
      numericInput(ns("min_count"), "Pre-filter: minimum count", value = 10, min = 0, step = 1),
      numericInput(ns("min_samples"), "... in at least N samples", value = 3, min = 1, step = 1),
      selectInput(ns("transform"), "Transformation for QC / heatmaps",
                  choices = c("Auto (rlog if <= 30 samples, else VST)" = "auto",
                              "VST" = "vst", "rlog" = "rlog")),
      actionButton(ns("run"), tagList(icon("play"), "Run DESeq2"),
                   class = "btn-primary", style = "width: 100%; margin-top: 10px;"),
      hr(),
      div(id = ns("log"), class = "deseq-log")
    ),
    mainPanel(
      width = 8,
      uiOutput(ns("status")),
      tabsetPanel(
        tabPanel("Counts", DT::dataTableOutput(ns("counts_table"))),
        tabPanel("Design", DT::dataTableOutput(ns("design_table")))
      )
    )
  )
}

mod_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    inputs <- reactiveVal(NULL)
    result <- reactiveVal(NULL)

    load_inputs <- function(counts_path, design_path, label) {
      tryCatch({
        cr <- read_counts(counts_path)
        dr <- read_design(design_path)
        counts <- match_samples(cr$counts, dr$design)
        if (length(usable_variables(dr$design)) == 0) {
          stop("No column of the design matrix has at least two levels.")
        }
        inputs(list(counts = counts, design = dr$design, label = label))
        result(NULL)
        for (m in c(cr$messages, dr$messages)) showNotification(m, type = "warning", duration = 15)
        showNotification(sprintf("Loaded %s: %d genes x %d samples.", label,
                                 nrow(counts), ncol(counts)), type = "message")
      }, error = function(e) {
        inputs(NULL)
        result(NULL)
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      })
    }

    observeEvent(list(input$counts_file, input$design_file), {
      req(input$counts_file, input$design_file)
      load_inputs(input$counts_file$datapath, input$design_file$datapath, "uploaded files")
    })

    observeEvent(input$load_example, {
      load_inputs(EXAMPLE_COUNTS, EXAMPLE_DESIGN, "example dataset (pasilla)")
    })

    output$model_ui <- renderUI({
      d <- inputs()
      if (is.null(d)) return(helpText("Load data to choose the model."))
      vars <- usable_variables(d$design)
      tagList(
        selectInput(ns("variable"), "Variable of interest", choices = vars, selected = vars[1]),
        selectizeInput(ns("covariates"), "Covariates (batch, replicate, ...)",
                       choices = vars, multiple = TRUE)
      )
    })

    # Keep covariates consistent with the variable of interest and propose a
    # sensible pre-filter (smallest group size).
    observeEvent(list(inputs(), input$variable), {
      d <- inputs()
      req(d, input$variable %in% names(d$design))
      vars <- setdiff(usable_variables(d$design), input$variable)
      updateSelectizeInput(session, "covariates", choices = vars,
                           selected = intersect(input$covariates, vars))
      updateNumericInput(session, "min_samples",
                         value = min(table(d$design[[input$variable]])),
                         max = ncol(d$counts))
    })

    output$status <- renderUI({
      d <- inputs()
      r <- result()
      if (is.null(d)) {
        return(div(class = "alert alert-info",
                   "Upload a count matrix and a design matrix, or load the example dataset."))
      }
      if (is.null(r)) {
        return(div(class = "alert alert-warning",
                   sprintf("%s loaded (%d genes x %d samples). Choose the model and click 'Run DESeq2'.",
                           d$label, nrow(d$counts), ncol(d$counts))))
      }
      pf <- S4Vectors::metadata(r$dds)$prefilter
      div(class = "alert alert-success",
          HTML(sprintf("DESeq2 done. Design: <code>%s</code>. %d / %d genes kept after pre-filtering. Transformation: %s.",
                       htmltools::htmlEscape(deparse(DESeq2::design(r$dds))),
                       pf$n_after, pf$n_before, S4Vectors::metadata(r$transformed)$method)))
    })

    output$counts_table <- DT::renderDataTable({
      d <- inputs()
      req(d)
      DT::datatable(as.data.frame(d$counts), options = list(pageLength = 10, scrollX = TRUE))
    })

    output$design_table <- DT::renderDataTable({
      d <- inputs()
      req(d)
      DT::datatable(d$design, options = list(pageLength = 10, scrollX = TRUE))
    })

    observeEvent(input$run, {
      d <- inputs()
      req(d, input$variable)
      shinyjs::html("log", "")
      log_message <- function(m) {
        shinyjs::html("log", paste0(htmltools::htmlEscape(conditionMessage(m)), "<br>"), add = TRUE)
        invokeRestart("muffleMessage")
      }
      withProgress(message = "Running DESeq2", value = 0, {
        tryCatch({
          incProgress(0.1, detail = "Fitting the model...")
          dds <- withCallingHandlers(
            run_deseq(d$counts, d$design, input$variable, input$covariates,
                      min_count = input$min_count, min_samples = input$min_samples),
            message = log_message)
          incProgress(0.6, detail = "Transforming counts...")
          transformed <- withCallingHandlers(transform_counts(dds, input$transform),
                                             message = log_message)
          result(list(dds = dds, transformed = transformed, design = d$design,
                      variable = input$variable, covariates = input$covariates))
          showNotification("DESeq2 analysis completed.", type = "message")
        }, error = function(e) {
          result(NULL)
          showNotification(paste("DESeq2 failed:", conditionMessage(e)),
                           type = "error", duration = NULL)
        })
      })
    })

    result
  })
}
