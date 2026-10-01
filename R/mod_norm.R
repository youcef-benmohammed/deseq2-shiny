# Module: normalised counts table and download.

mod_norm_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText("Counts divided by DESeq2 size factors (median-of-ratios).",
               "Use them for visualisation, not as input to another DE tool."),
      downloadButton(ns("download"), "Download normalised counts", class = "btn-success"),
      hr(),
      tableOutput(ns("size_factors"))
    ),
    mainPanel(width = 9, DT::dataTableOutput(ns("table")))
  )
}

mod_norm_server <- function(id, result) {
  moduleServer(id, function(input, output, session) {
    normalized <- reactive({
      r <- result()
      req(r)
      DESeq2::counts(r$dds, normalized = TRUE)
    })

    output$table <- DT::renderDataTable({
      DT::datatable(round(normalized(), 2), options = list(pageLength = 10, scrollX = TRUE))
    })

    output$size_factors <- renderTable({
      r <- result()
      req(r)
      data.frame(Sample = colnames(r$dds), `Size factor` = round(DESeq2::sizeFactors(r$dds), 3),
                 check.names = FALSE)
    })

    output$download <- downloadHandler(
      filename = function() "normalized_counts.tsv",
      content = function(file) {
        mat <- normalized()
        write_tsv(cbind(gene = rownames(mat), as.data.frame(round(mat, 3))), file)
      }
    )
  })
}
