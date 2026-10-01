# Module: differential expression (contrast, plots, tables, report).

STATUS_COLORS <- c("Up" = "#C0392B", "Down" = "#2471A3",
                   "Not significant" = "#B3B3B3", "Not tested (padj = NA)" = "#E0E0E0")

mod_de_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      uiOutput(ns("contrast_ui")),
      hr(),
      numericInput(ns("alpha"), "Adjusted p-value (FDR) threshold", value = 0.05,
                   min = 0.001, max = 0.2, step = 0.01),
      numericInput(ns("lfc"), "|log2 fold change| threshold", value = 1, min = 0, max = 10, step = 0.25),
      checkboxInput(ns("shrink"), "Shrink log2 fold changes (lfcShrink, type 'normal')", value = TRUE),
      helpText("Shrinkage reduces the noisy fold changes of low-count genes;",
               "the unshrunk estimate is kept in the 'log2FoldChange_MLE' column."),
      actionButton(ns("run"), tagList(icon("magic"), "Run DE analysis"),
                   class = "btn-primary", style = "width: 100%;"),
      hr(),
      numericInput(ns("n_top"), "Heatmap: number of top DE genes", value = 50, min = 2, max = 500),
      hr(),
      downloadButton(ns("download_sig"), "Significant genes (TSV)"),
      br(), br(),
      downloadButton(ns("download_all"), "All genes (TSV)"),
      br(), br(),
      downloadButton(ns("report"), "HTML report", class = "btn-success")
    ),
    mainPanel(
      width = 9,
      uiOutput(ns("placeholder")),
      fluidRow(
        column(5, h4("Summary"), tableOutput(ns("summary"))),
        column(7, h4("Gene detail"), helpText("Click a point in the volcano / MA plot or a table row."),
               plotOutput(ns("gene_plot"), height = "260px"))
      ),
      hr(),
      fluidRow(
        column(6, plotlyOutput(ns("volcano"), height = "450px")),
        column(6, plotlyOutput(ns("ma"), height = "450px"))
      ),
      hr(),
      fluidRow(
        column(6, plotlyOutput(ns("pvalues"), height = "320px")),
        column(6, helpText("A well-behaved test shows a flat distribution with a peak near 0.",
                           "A peak near 1 or a U-shape suggests a mis-specified model or hidden batch effects."))
      ),
      hr(),
      plotlyOutput(ns("heatmap"), height = "650px"),
      hr(),
      h4("Significant genes"),
      DT::dataTableOutput(ns("table"))
    )
  )
}

mod_de_server <- function(id, result) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    res_raw <- reactiveVal(NULL)
    selected_gene <- reactiveVal(NULL)

    # A new DESeq2 run invalidates previous DE results.
    observeEvent(result(), {
      res_raw(NULL)
      selected_gene(NULL)
    }, ignoreNULL = FALSE)

    output$placeholder <- renderUI({
      if (is.null(result())) {
        div(class = "alert alert-info", "Run DESeq2 first (tab 'Data & model').")
      } else if (is.null(res_raw())) {
        div(class = "alert alert-info", "Choose the contrast and click 'Run DE analysis'.")
      }
    })

    output$contrast_ui <- renderUI({
      r <- result()
      req(r)
      lv <- levels(r$dds[[r$variable]])
      tagList(
        p(strong("Variable:"), code(r$variable)),
        selectInput(ns("numerator"), "Numerator (e.g. treated)", choices = lv, selected = lv[2]),
        selectInput(ns("denominator"), "Denominator / reference (e.g. control)",
                    choices = lv, selected = lv[1])
      )
    })

    observeEvent(input$run, {
      r <- result()
      req(r, input$numerator, input$denominator)
      if (input$numerator == input$denominator) {
        showNotification("Numerator and denominator must be different.", type = "error")
        return()
      }
      withProgress(message = "Computing DE results...", value = 0.3, {
        tryCatch({
          df <- get_results(r$dds, r$variable, input$numerator, input$denominator,
                            alpha = input$alpha, lfc_threshold = input$lfc,
                            shrink = input$shrink)
          res_raw(list(df = df, numerator = input$numerator, denominator = input$denominator,
                       shrink = input$shrink, alpha_filtering = input$alpha))
          selected_gene(df$gene[1])
        }, error = function(e) {
          showNotification(paste("DE analysis failed:", conditionMessage(e)),
                           type = "error", duration = NULL)
        })
      })
    })

    # Thresholds can be changed without re-running DESeq2.
    res <- reactive({
      x <- res_raw()
      req(x, input$alpha, !is.null(input$lfc))
      x$df$status <- classify_genes(x$df$padj, x$df$log2FoldChange, input$alpha, input$lfc)
      x
    })

    contrast_label <- reactive(sprintf("%s vs %s", res()$numerator, res()$denominator))

    output$summary <- renderTable(summarise_results(res()$df, result()$dds), digits = 0)

    output$volcano <- renderPlotly({
      x <- res()
      df <- x$df[!is.na(x$df$padj), ]
      df$neglog10 <- -log10(pmax(df$padj, .Machine$double.xmin))
      plot_ly(df, x = ~log2FoldChange, y = ~neglog10, color = ~status, colors = STATUS_COLORS,
              key = ~gene, text = ~gene, type = "scattergl", mode = "markers",
              marker = list(size = 5, opacity = 0.7), source = ns("plots"),
              hovertemplate = "%{text}<br>log2FC: %{x:.2f}<br>-log10 padj: %{y:.2f}<extra></extra>") %>%
        layout(title = paste("Volcano plot:", contrast_label()),
               xaxis = list(title = "log2 fold change"),
               yaxis = list(title = "-log10 adjusted p-value"),
               shapes = list(
                 vline(-input$lfc), vline(input$lfc),
                 hline(-log10(input$alpha)))) %>%
        event_register("plotly_click")
    })

    output$ma <- renderPlotly({
      x <- res()
      df <- x$df[x$df$baseMean > 0, ]
      plot_ly(df, x = ~baseMean, y = ~log2FoldChange, color = ~status, colors = STATUS_COLORS,
              key = ~gene, text = ~gene, type = "scattergl", mode = "markers",
              marker = list(size = 5, opacity = 0.7), source = ns("plots"),
              hovertemplate = "%{text}<br>mean: %{x:.1f}<br>log2FC: %{y:.2f}<extra></extra>") %>%
        layout(title = paste("MA plot:", contrast_label()),
               xaxis = list(title = "Mean of normalised counts", type = "log"),
               yaxis = list(title = "log2 fold change"),
               shapes = list(hline(-input$lfc), hline(input$lfc), hline(0, dash = "solid"))) %>%
        event_register("plotly_click")
    })

    observe({
      req(res_raw())
      # Plots may not be rendered yet on the first run; silence plotly's warning.
      ev <- suppressWarnings(event_data("plotly_click", source = ns("plots")))
      if (!is.null(ev$key)) selected_gene(as.character(ev$key[[1]]))
    })

    output$pvalues <- renderPlotly({
      x <- res()
      df <- x$df[!is.na(x$df$pvalue), ]
      plot_ly(df, x = ~pvalue, type = "histogram", xbins = list(start = 0, end = 1, size = 0.05),
              marker = list(color = "#2c3e50", line = list(color = "white", width = 1))) %>%
        layout(title = "Raw p-value distribution", xaxis = list(title = "p-value"),
               yaxis = list(title = "Genes"))
    })

    output$heatmap <- renderPlotly({
      x <- res()
      r <- result()
      z <- top_gene_matrix(r$transformed, x$df, n = input$n_top)
      validate(need(!is.null(z), "Fewer than 2 significant genes: nothing to show in the heatmap."))
      heatmaply::heatmaply(z, col_side_colors = r$design[colnames(z), , drop = FALSE],
                           scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                             low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0),
                           main = sprintf("Top %d DE genes (row z-score, %s)", nrow(z),
                                          S4Vectors::metadata(r$transformed)$method),
                           showticklabels = c(TRUE, nrow(z) <= 30))
    })

    significant <- reactive({
      df <- res()$df
      df[df$status %in% c("Up", "Down"), ]
    })

    output$table <- DT::renderDataTable({
      df <- significant()
      num <- vapply(df, is.numeric, logical(1))
      DT::datatable(df, selection = "single", rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) %>%
        DT::formatSignif(names(df)[num], digits = 4)
    })

    observeEvent(input$table_rows_selected, {
      selected_gene(significant()$gene[input$table_rows_selected])
    })

    output$gene_plot <- renderPlot({
      r <- result()
      g <- selected_gene()
      req(r, g, g %in% rownames(r$dds))
      d <- DESeq2::plotCounts(r$dds, gene = g, intgroup = r$variable, returnData = TRUE)
      ggplot2::ggplot(d, ggplot2::aes(x = .data[[r$variable]], y = count, colour = .data[[r$variable]])) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, height = 0), size = 3) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(title = g, x = NULL, y = "Normalised count (log10)") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "none")
    })

    file_stub <- function() {
      make.names(sprintf("DESeq2_%s_vs_%s", res()$numerator, res()$denominator))
    }

    output$download_sig <- downloadHandler(
      filename = function() paste0(file_stub(), "_significant.tsv"),
      content = function(file) write_tsv(significant(), file)
    )

    output$download_all <- downloadHandler(
      filename = function() paste0(file_stub(), "_all_genes.tsv"),
      content = function(file) write_tsv(res()$df, file)
    )

    output$report <- downloadHandler(
      filename = function() paste0(file_stub(), "_report.html"),
      content = function(file) {
        withProgress(message = "Rendering report...", value = 0.5, {
          tmp <- tempfile(fileext = ".Rmd")
          file.copy(REPORT_TEMPLATE, tmp, overwrite = TRUE)
          r <- result()
          x <- res()
          rmarkdown::render(
            tmp, output_file = file, quiet = TRUE, envir = new.env(parent = environment()),
            params = list(
              dds = r$dds, transformed = r$transformed, design = r$design,
              variable = r$variable, results = x$df,
              numerator = x$numerator, denominator = x$denominator,
              alpha = input$alpha, lfc = input$lfc, shrink = x$shrink,
              n_top = input$n_top))
        })
      }
    )
  })
}

vline <- function(x) {
  list(type = "line", x0 = x, x1 = x, yref = "paper", y0 = 0, y1 = 1,
       line = list(dash = "dash", color = "#555555", width = 1))
}

hline <- function(y, dash = "dash") {
  list(type = "line", y0 = y, y1 = y, xref = "paper", x0 = 0, x1 = 1,
       line = list(dash = dash, color = "#555555", width = 1))
}
