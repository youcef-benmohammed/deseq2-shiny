# Run from the repository root:  Rscript tests/testthat.R
library(testthat)
suppressPackageStartupMessages(library(DESeq2))
source("R/helpers.R")
test_dir("tests/testthat", stop_on_failure = TRUE)
