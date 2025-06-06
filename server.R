library(shiny)
library(DT)
library(shinyjs)
library(plotly)
library(DESeq2)
library(heatmaply)
library(ggplot2)
library(RColorBrewer)
shinyServer(function(input, output, session) {
  # Reactive values to store DESeqDataSet and normalized counts
  dds_r <- reactiveVal(NULL)
  normalized_counts_r <- reactiveVal(NULL)
  
  # Load counts data
  counts_data <- reactive({
    req(input$counts_file)
    df <- read.table(
      input$counts_file$datapath,
      header = TRUE,
      sep = "\t",
      row.names = 1,
      check.names = FALSE
    )
    # Convert to integer matrix
    mat <- as.matrix(df)
    storage.mode(mat) <- "integer"
    return(mat)
  })
  
  # Load design matrix
  design_data <- reactive({
    req(input$design_matrix)
    df <- read.table(
      input$design_matrix$datapath,
      header = TRUE,
      sep = "\t",
      row.names = 1,
      check.names = FALSE
    )
    # Use sample column as rownames if present, then remove it
    if ("sample" %in% colnames(df)) {
      rownames(df) <- df$sample
      df$sample <- NULL
    }
    # Convert all columns to factors
    df[] <- lapply(df, as.factor)
    return(df)
  })
  
  # Render uploaded tables
  output$counts_table <- DT::renderDataTable({
    DT::datatable(counts_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$design_matrix_table <- DT::renderDataTable({
    DT::datatable(design_data(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  # Run DESeq analysis when button clicked
  observeEvent(input$run, {
    req(counts_data(), design_data())
    
    # Check sample names match
    if (!all(colnames(counts_data()) %in% rownames(design_data()))) {
      showNotification("Sample names don't match.", type = "error")
      return(NULL)
    }
    
    progress <- Progress$new(session, min = 1, max = 100)
    on.exit(progress$close())
    progress$set(message = 'DESeq analysis...', value = 0)
    
    counts <- counts_data()
    design <- design_data()
    counts <- counts[, rownames(design)]
    progress$set(detail = 'Building DESeqDataSet...', value = 20)
    
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = design,
      design = ~ genotype + replicate
    )
    
    progress$set(detail = "Running DESeq()...", value = 50)
    
    # Run DESeq with message handlers to display logs if needed
    withCallingHandlers({
      dds_res <- DESeq(dds)
      dds_r(dds_res)
    }, message = function(m) {
      shinyjs::html("log", paste0(m$message, "<br>"), add = TRUE)
    })
    
    # Calculate normalized counts
    dds_r_current <- dds_r()
    normalized_counts_r(counts(dds_r_current, normalized = TRUE))
    
    progress$set(detail = "Preparing QC plots and UI...", value = 80)
    
    # Prepare UI elements dynamically
    group_cols <- colnames(colData(dds_r_current))
    
    output$filter_col_select <- renderUI({
      selectInput("filter_col",
                  "Filter by:",
                  choices = group_cols,
                  selected = group_cols[1])
    })
    
    output$filter_values_ui <- renderUI({
      req(input$filter_col)
      vals <- unique(as.character(colData(dds_r_current)[[input$filter_col]]))
      selectInput(
        "filter_values",
        paste0("Select values for ", input$filter_col),
        choices = vals,
        selected = vals,
        multiple = TRUE,
        selectize = TRUE,
        width = "100%"
      )
    })
    
    output$group_by_select <- renderUI({
      selectInput("group_by", "Group by", choices = group_cols)
    })
    
    output$pc_x_select <- renderUI({
      selectInput(
        "pc_x",
        "Select x-axis PC",
        choices = paste0("PC", 1:10),
        # max 10 PCs shown by default
        selected = "PC1"
      )
    })
    
    output$pc_y_select <- renderUI({
      selectInput(
        "pc_y",
        "Select y-axis PC",
        choices = paste0("PC", 1:10),
        selected = "PC2"
      )
    })
    
    # Compute rlog and PCA for plots
    rld <- rlog(dds_r_current, blind = TRUE)
    rld_mat <- assay(rld)
    pca_res <- prcomp(t(rld_mat), center = TRUE, scale. = FALSE)
    percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
    
    pca_data <- data.frame(pca_res$x)
    pca_data$name <- colnames(rld_mat)
    pca_data <- cbind(pca_data, as.data.frame(colData(dds_r_current)))
    
    # Store PCA data in reactiveVal to avoid recomputation in plots
    pca_data_r <- reactiveVal(pca_data)
    percentVar_r <- reactiveVal(percentVar)
    
    # PCA plot
    output$pcaPlot <- renderPlotly({
      req(
        input$filter_col,
        input$filter_values,
        input$group_by,
        input$pc_x,
        input$pc_y
      )
      filtered_data <- pca_data_r()[pca_data_r()[[input$filter_col]] %in% input$filter_values, ]
      x_pc <- input$pc_x
      y_pc <- input$pc_y
      x_idx <- as.numeric(sub("PC", "", x_pc))
      y_idx <- as.numeric(sub("PC", "", y_pc))
      
      plot_ly(
        data = filtered_data,
        x = ~ get(x_pc),
        y = ~ get(y_pc),
        type = 'scatter',
        mode = 'markers',
        color = ~ get(input$group_by),
        text = ~ name,
        marker = list(size = 10)
      ) %>%
        layout(
          title = "PCA plot",
          xaxis = list(
            title = paste0(x_pc, ": ", percentVar_r()[x_idx], "% variance"),
            zeroline = FALSE
          ),
          yaxis = list(
            title = paste0(y_pc, ": ", percentVar_r()[y_idx], "% variance"),
            zeroline = FALSE
          )
        )
    })
    
    # Elbow plot for variance explained
    output$elbowPlot <- renderPlotly({
      var_df <- data.frame(PC = seq_along(percentVar_r()), Variance = percentVar_r())
      plot_ly(
        data = var_df,
        x = ~ PC,
        y = ~ Variance,
        type = 'scatter',
        mode = 'lines+markers',
        marker = list(size = 8, color = 'blue'),
        line = list(color = 'blue')
      ) %>%
        layout(
          title = "Elbow plot",
          xaxis = list(title = "PC"),
          yaxis = list(title = "Percent Variance", range = c(0, 100))
        )
    })
    
    # Sample clustering heatmap
    rld_cor <- cor(assay(rld))
    output$sampleClustering <- renderPlotly({
      heatmaply(rld_cor, main = "Hierarchical Clustering of Samples (rlog)")
    })
    
    # Normalized counts table and download
    output$normalized_counts <- DT::renderDataTable({
      DT::datatable(normalized_counts_r(),
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$normcounts <- downloadHandler(
      filename = function() {
        "normalizedCounts.tsv"
      },
      content = function(file) {
        write.table(
          normalized_counts_r(),
          file,
          sep = "\t",
          quote = FALSE,
          row.names = TRUE,
          col.names = NA
        )
      }
    )
    
  })
  
  
  # DEGs analysis reactive
  dds.results <- eventReactive(input$run_deg, {
    req(dds_r(), input$cond_a, input$cond_b)
    results(
      dds_r(),
      contrast = c("genotype", input$cond_b, input$cond_a),
      pAdjustMethod = "fdr"
    )
  })
  
  # UI for conditions
  output$condition_a <- renderUI({
    req(dds_r())
    selectInput("cond_a", "Condition A", choices = levels(colData(dds_r())$genotype))
  })
  
  output$condition_b <- renderUI({
    req(dds_r(), input$cond_a)
    choices <- setdiff(levels(colData(dds_r())$genotype), input$cond_a)
    selectInput("cond_b", "Condition B", choices = choices)
  })
  
  # summary
  output$degSummary <- renderPrint({
    req(dds.results())
    summary(dds.results())
  })
  
  # --------------------------- MA plot plot ----------------------------
  output$MAPlot <- renderPlotly({
    req(dds.results(), input$logfc_thresh, input$padj_thresh)
    
    res <- dds.results()
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    df <- na.omit(df)
    df$mean <- df$baseMean + 1
    df$significant <- with(df,
                           padj < input$padj_thresh &
                             abs(log2FoldChange) > input$logfc_thresh)
    
    colors <- c("grey", "red")
    
    p <- ggplot(df,
                aes(
                  x = log10(mean),
                  y = log2FoldChange,
                  color = significant,
                  text = gene
                )) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = colors) +
      geom_hline(
        yintercept = c(-input$logfc_thresh, input$logfc_thresh),
        linetype = "dashed",
        color = "blue"
      ) +
      geom_hline(yintercept = 0, linetype = "solid") +
      labs(title = "MA Plot", x = "Log10 mean normalized count", y = "Log2 Fold Change") +
      theme_minimal()
    
    ggplotly(p, tooltip = c("text", "x", "y"))
  })
  
  # --------------------------- P-values plot ----------------------------
  output$pvalues <- renderPlotly({
    req(dds.results())
    res <- as.data.frame(dds.results())
    res <- res[!is.na(res$pvalue), ]
    
    p <- ggplot(res, aes(x = pvalue)) +
      geom_histogram(bins = 20, fill = "#2c3e50", color = "white") +
      labs(title = "p-values distribution",
           x = "p-value",
           y = "count") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # --------------------------- Volcano plot ----------------------------
  output$volcanoPlot <- renderPlotly({
    req(dds.results(), input$logfc_thresh, input$padj_thresh)
    
    res <- dds.results()
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    df <- na.omit(df)
    
    df$log10padj <- -log10(df$padj)
    df$log10padj[is.infinite(df$log10padj)] <- max(df$log10padj[is.finite(df$log10padj)], na.rm =
                                                     TRUE) + 1
    
    df$significant <- with(df,
                           padj < input$padj_thresh &
                             abs(log2FoldChange) > input$logfc_thresh)
    
    colors <- c("grey", "red")
    
    p <- ggplot(df,
                aes(
                  x = log2FoldChange,
                  y = log10padj,
                  color = significant,
                  text = gene
                )) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = colors) +
      geom_vline(
        xintercept = c(-input$logfc_thresh, input$logfc_thresh),
        linetype = "dashed",
        color = "blue"
      ) +
      geom_hline(
        yintercept = -log10(input$padj_thresh),
        linetype = "dashed",
        color = "blue"
      ) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
      theme_minimal()
    
    ggplotly(p, tooltip = c("text", "x", "y"))
  })

  # --------------------------- Table of DEGs ----------------------------
  output$degTable <- DT::renderDataTable({
    req(dds.results(), input$logfc_thresh, input$padj_thresh)
    
    df <- as.data.frame(dds.results())
    df$gene <- rownames(df)
    df <- na.omit(df)
    
    filtered <- subset(df,
                       padj < input$padj_thresh &
                         abs(log2FoldChange) > input$logfc_thresh)
    
    DT::datatable(filtered, options = list(pageLength = 10, scrollX = TRUE))
  })
  # --------------------------- Download Table of DEGs ----------------------------
  output$downloadDEG <- downloadHandler(
    filename = function() {
      paste0("DEGs_", input$cond_b, "_vs_", input$cond_a, ".tsv")
    },
    content = function(file) {
      df <- as.data.frame(dds.results())
      df$gene <- rownames(df)
      df <- na.omit(df)
      filtered <- subset(df,
                         padj < input$padj_thresh &
                           abs(log2FoldChange) > input$logfc_thresh)
      write.table(
        filtered,
        file,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE
      )
    }
  )
  
})

