# Loaded once at app start-up. Files in R/ are sourced automatically by Shiny.
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinythemes)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(heatmaply)
  library(DESeq2)
})

options(shiny.maxRequestSize = 200 * 1024^2)  # allow uploads up to 200 MB

EXAMPLE_COUNTS <- "inst/extdata/pasilla_counts.tsv"
EXAMPLE_DESIGN <- "inst/extdata/pasilla_design.tsv"
REPORT_TEMPLATE <- normalizePath("report/report.Rmd")
