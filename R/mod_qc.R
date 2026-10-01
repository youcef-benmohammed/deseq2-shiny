# Module: quality-control plots (PCA, scree plot, sample correlation).

mod_qc_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      uiOutput(ns("controls")),
      numericInput(ns("ntop"), "PCA: number of most variable genes", value = 500, min = 50, step = 50)
    ),
    mainPanel(
      width = 9,
      uiOutput(ns("placeholder")),
      plotlyOutput(ns("pca"), height = "450px"),
      hr(),
      plotlyOutput(ns("scree"), height = "300px"),
      hr(),
      plotlyOutput(ns("correlation"), height = "600px")
    )
  )
}

mod_qc_server <- function(id, result) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$placeholder <- renderUI({
      if (is.null(result())) div(class = "alert alert-info", "Run DESeq2 first (tab 'Data & model').")
    })

    pca <- reactive({
      r <- result()
      req(r, input$ntop)
      compute_pca(SummarizedExperiment::assay(r$transformed), r$design, ntop = input$ntop)
    })

    output$controls <- renderUI({
      r <- result()
      req(r)
      vars <- names(r$design)
      pcs <- paste0("PC", seq_len(ncol(r$dds)))
      pcs <- pcs[seq_len(min(length(pcs), 10))]
      tagList(
        selectInput(ns("color_by"), "Colour by", choices = vars, selected = r$variable),
        selectInput(ns("shape_by"), "Shape by", choices = c("(none)", vars)),
        selectInput(ns("filter_col"), "Show samples where", choices = vars, selected = r$variable),
        uiOutput(ns("filter_values_ui")),
        fluidRow(
          column(6, selectInput(ns("pc_x"), "X axis", choices = pcs, selected = pcs[1])),
          column(6, selectInput(ns("pc_y"), "Y axis", choices = pcs, selected = pcs[min(2, length(pcs))]))
        )
      )
    })

    output$filter_values_ui <- renderUI({
      r <- result()
      req(r, input$filter_col %in% names(r$design))
      vals <- levels(r$design[[input$filter_col]])
      selectizeInput(ns("filter_values"), NULL, choices = vals, selected = vals, multiple = TRUE)
    })

    output$pca <- renderPlotly({
      p <- pca()
      req(input$color_by %in% names(p$data), input$pc_x %in% names(p$data),
          input$pc_y %in% names(p$data), input$filter_col %in% names(p$data))
      df <- p$data[p$data[[input$filter_col]] %in% input$filter_values, , drop = FALSE]
      validate(need(nrow(df) > 0, "No sample selected."))
      ix <- as.integer(sub("PC", "", input$pc_x))
      iy <- as.integer(sub("PC", "", input$pc_y))
      args <- list(data = df, x = df[[input$pc_x]], y = df[[input$pc_y]],
                   type = "scatter", mode = "markers",
                   color = df[[input$color_by]], colors = palette_for(df[[input$color_by]]),
                   text = df$sample,
                   hoverinfo = "text", marker = list(size = 11))
      if (!identical(input$shape_by, "(none)") && input$shape_by %in% names(df)) {
        args$symbol <- df[[input$shape_by]]
      }
      do.call(plot_ly, args) %>%
        layout(title = "PCA",
               xaxis = list(title = sprintf("%s: %s%% variance", input$pc_x, p$percent[ix]), zeroline = FALSE),
               yaxis = list(title = sprintf("%s: %s%% variance", input$pc_y, p$percent[iy]), zeroline = FALSE))
    })

    output$scree <- renderPlotly({
      p <- pca()
      df <- data.frame(PC = seq_along(p$percent), Variance = p$percent)
      plot_ly(df, x = ~PC, y = ~Variance, type = "scatter", mode = "lines+markers") %>%
        layout(title = "Variance explained per component (scree plot)",
               xaxis = list(title = "Principal component", dtick = 1),
               yaxis = list(title = "% variance", range = c(0, 100)))
    })

    output$correlation <- renderPlotly({
      r <- result()
      req(r)
      mat <- stats::cor(SummarizedExperiment::assay(r$transformed))
      heatmaply::heatmaply(mat, row_side_colors = r$design[colnames(mat), , drop = FALSE],
                           main = sprintf("Sample correlation (%s)", S4Vectors::metadata(r$transformed)$method),
                           limits = c(min(mat), 1))
    })
  })
}
