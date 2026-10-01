# Pure helper functions (no Shiny dependency) so they can be unit-tested.

# ---------------------------------------------------------------------------
# Input reading and validation
# ---------------------------------------------------------------------------

#' Read a tab-delimited count matrix (genes in rows, samples in columns).
#'
#' The first column must contain gene IDs. Non-integer values (e.g. Salmon /
#' RSEM estimated counts) are rounded and a warning message is returned.
#'
#' @return list(counts = integer matrix, messages = character vector)
read_counts <- function(path) {
  df <- utils::read.delim(path, row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE)
  if (ncol(df) < 2) {
    stop("The count matrix must contain at least two sample columns ",
         "(is the file tab-delimited?).")
  }
  non_numeric <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("Non-numeric column(s) in the count matrix: ",
         paste(utils::head(non_numeric, 5), collapse = ", "))
  }
  mat <- as.matrix(df)
  if (anyNA(mat)) stop("The count matrix contains missing values (NA).")
  if (any(mat < 0)) stop("The count matrix contains negative values.")

  messages <- character(0)
  if (any(mat != round(mat))) {
    mat <- round(mat)
    messages <- c(messages, paste(
      "Non-integer counts detected (e.g. Salmon/RSEM estimates):",
      "values were rounded to the nearest integer. For transcript-level",
      "quantifications, prefer tximport + DESeqDataSetFromTximport."))
  }
  if (any(mat > .Machine$integer.max)) {
    stop("Some counts exceed the maximum integer value supported by R.")
  }
  storage.mode(mat) <- "integer"
  list(counts = mat, messages = messages)
}

#' Read a tab-delimited sample sheet. The first column must contain sample
#' IDs; every other column is converted to a factor.
#'
#' @return list(design = data.frame, messages = character vector)
read_design <- function(path) {
  df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = "character")
  if (ncol(df) < 2) {
    stop("The design matrix must contain a sample ID column and at least ",
         "one variable column (is the file tab-delimited?).")
  }
  ids <- df[[1]]
  if (anyDuplicated(ids)) {
    stop("Duplicated sample IDs in the design matrix: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))
  }
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids

  messages <- character(0)
  clean <- make.names(colnames(df), unique = TRUE)
  if (!identical(clean, colnames(df))) {
    messages <- c(messages, paste0(
      "Column names were made syntactically valid for the design formula: ",
      paste(sprintf("'%s' -> '%s'", colnames(df)[clean != colnames(df)],
                    clean[clean != colnames(df)]), collapse = ", ")))
    colnames(df) <- clean
  }
  df[] <- lapply(df, function(x) factor(x, levels = unique(x)))
  list(design = df, messages = messages)
}

#' Check that sample IDs agree between counts and design, and reorder the
#' count columns to follow the design rows.
match_samples <- function(counts, design) {
  only_counts <- setdiff(colnames(counts), rownames(design))
  only_design <- setdiff(rownames(design), colnames(counts))
  if (length(only_counts) || length(only_design)) {
    fmt <- function(x) paste(utils::head(x, 10), collapse = ", ")
    msg <- "Sample IDs do not match between the count and design matrices."
    if (length(only_counts)) msg <- paste0(msg, "\n- Only in counts: ", fmt(only_counts))
    if (length(only_design)) msg <- paste0(msg, "\n- Only in design: ", fmt(only_design))
    stop(msg)
  }
  counts[, rownames(design), drop = FALSE]
}

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

#' Variables of the design that can be used in a model (>= 2 levels).
usable_variables <- function(design) {
  names(design)[vapply(design, function(x) nlevels(droplevels(x)) >= 2, logical(1))]
}

#' Build a DESeq2 design formula. The variable of interest is placed last,
#' as recommended by DESeq2.
build_formula <- function(variable, covariates = character(0)) {
  covariates <- setdiff(covariates, variable)
  stats::as.formula(paste("~", paste(c(covariates, variable), collapse = " + ")))
}

#' Run the full DESeq2 workflow.
#'
#' @param min_count,min_samples pre-filtering: keep genes with at least
#'   `min_count` reads in at least `min_samples` samples.
run_deseq <- function(counts, design, variable, covariates = character(0),
                      min_count = 10, min_samples = NULL) {
  if (!variable %in% names(design)) stop("Unknown variable: ", variable)
  bad <- setdiff(c(variable, covariates), usable_variables(design))
  if (length(bad)) stop("These variables have a single level: ", paste(bad, collapse = ", "))

  if (is.null(min_samples)) min_samples <- min(table(design[[variable]]))
  formula <- build_formula(variable, covariates)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = design,
                                        design = formula)
  n_before <- nrow(dds)
  keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
  dds <- dds[keep, ]
  if (nrow(dds) == 0) stop("No gene passes the pre-filtering thresholds.")
  dds <- DESeq2::DESeq(dds)
  S4Vectors::metadata(dds)$prefilter <- list(
    min_count = min_count, min_samples = min_samples,
    n_before = n_before, n_after = nrow(dds))
  dds
}

#' Variance-stabilising transformation. "auto" uses rlog for <= 30 samples
#' and VST otherwise (rlog becomes very slow on large designs).
transform_counts <- function(dds, method = c("auto", "vst", "rlog")) {
  method <- match.arg(method)
  if (method == "auto") method <- if (ncol(dds) > 30) "vst" else "rlog"
  out <- if (method == "rlog") {
    DESeq2::rlog(dds, blind = TRUE)
  } else if (nrow(dds) < 1000) {
    DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  } else {
    DESeq2::vst(dds, blind = TRUE)
  }
  S4Vectors::metadata(out)$method <- method
  out
}

#' PCA on the `ntop` most variable genes (same approach as DESeq2::plotPCA).
compute_pca <- function(mat, coldata, ntop = 500) {
  vars <- matrixStats::rowVars(mat)
  select <- order(vars, decreasing = TRUE)[seq_len(min(ntop, length(vars)))]
  pca <- stats::prcomp(t(mat[select, , drop = FALSE]), center = TRUE, scale. = FALSE)
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  data <- cbind(as.data.frame(pca$x), sample = colnames(mat), as.data.frame(coldata))
  list(data = data, percent = percent, n_pc = ncol(pca$x))
}

#' Differential expression results for `numerator` vs `denominator`.
#'
#' Returns all genes (including those with padj = NA, which DESeq2 sets for
#' outliers and independently-filtered low counts) with a `status` column.
#' When `shrink = TRUE`, `log2FoldChange` is the shrunken estimate
#' (lfcShrink, type "normal") and the MLE is kept in `log2FoldChange_MLE`.
get_results <- function(dds, variable, numerator, denominator,
                        alpha = 0.05, lfc_threshold = 1, shrink = TRUE) {
  contrast <- c(variable, numerator, denominator)
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  df <- as.data.frame(res)
  if (shrink) {
    shr <- suppressMessages(
      DESeq2::lfcShrink(dds, contrast = contrast, res = res, type = "normal"))
    df$log2FoldChange_MLE <- df$log2FoldChange
    df$log2FoldChange <- shr$log2FoldChange
    df$lfcSE <- shr$lfcSE
  }
  df <- cbind(gene = rownames(df), df)
  rownames(df) <- NULL
  df$status <- classify_genes(df$padj, df$log2FoldChange, alpha, lfc_threshold)
  df[order(df$padj, na.last = TRUE), ]
}

classify_genes <- function(padj, lfc, alpha, lfc_threshold) {
  status <- ifelse(is.na(padj), "Not tested (padj = NA)",
             ifelse(padj < alpha & lfc > lfc_threshold, "Up",
               ifelse(padj < alpha & lfc < -lfc_threshold, "Down", "Not significant")))
  factor(status, levels = c("Up", "Down", "Not significant", "Not tested (padj = NA)"))
}

#' One-line-per-category summary used in the UI and the report.
summarise_results <- function(df, dds) {
  pf <- S4Vectors::metadata(dds)$prefilter
  tab <- table(df$status)
  data.frame(
    Category = c("Genes in input", "Genes after pre-filtering", names(tab)),
    Genes = c(pf$n_before, pf$n_after, as.integer(tab)),
    check.names = FALSE)
}

#' Write a data.frame as TSV with a correct header (no row names).
write_tsv <- function(df, file) {
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

#' Z-scored transformed expression of the `n` top DE genes.
top_gene_matrix <- function(transformed, res, n = 50) {
  sig <- res[res$status %in% c("Up", "Down"), ]
  genes <- utils::head(sig$gene, n)
  if (length(genes) < 2) return(NULL)
  mat <- SummarizedExperiment::assay(transformed)[genes, , drop = FALSE]
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

#' Qualitative palette with one colour per level (avoids RColorBrewer's
#' "minimal value for n is 3" warning with two groups).
palette_for <- function(x) {
  n <- nlevels(as.factor(x))
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  if (n <= length(pal)) pal[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}
