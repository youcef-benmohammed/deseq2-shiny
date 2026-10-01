run_example <- function(...) {
  counts <- read_counts(example_counts())$counts
  design <- read_design(example_design())$design
  counts <- match_samples(counts, design)
  suppressMessages(run_deseq(counts, design, ...))
}

test_that("full workflow runs on the example dataset with any column names", {
  dds <- run_example("condition", "type")
  expect_equal(deparse(design(dds)), "~type + condition")
  pf <- S4Vectors::metadata(dds)$prefilter
  expect_lt(pf$n_after, pf$n_before)

  res <- get_results(dds, "condition", "treated", "untreated", alpha = 0.05, lfc_threshold = 1)
  expect_equal(nrow(res), nrow(dds))
  expect_true(all(c("gene", "log2FoldChange", "log2FoldChange_MLE", "padj", "status") %in% names(res)))
  # pasilla knock-down: the pasilla gene itself must be strongly down-regulated
  expect_equal(as.character(res$status[res$gene == "FBgn0261552"]), "Down")
  expect_gt(sum(res$status %in% c("Up", "Down")), 50)

  s <- summarise_results(res, dds)
  expect_equal(s$Genes[1], pf$n_before)
  expect_equal(sum(s$Genes[-(1:2)]), pf$n_after)
})

test_that("run_deseq rejects single-level variables", {
  counts <- read_counts(example_counts())$counts
  design <- read_design(example_design())$design
  design$constant <- factor("x")
  expect_error(run_deseq(counts, design, "constant"), "single level")
})

test_that("transform_counts chooses rlog for small designs and VST for large ones", {
  small <- suppressMessages(DESeq(makeExampleDESeqDataSet(n = 1500, m = 6)))
  expect_equal(S4Vectors::metadata(transform_counts(small))$method, "rlog")
  large <- suppressMessages(DESeq(makeExampleDESeqDataSet(n = 1500, m = 32)))
  expect_equal(S4Vectors::metadata(transform_counts(large))$method, "vst")
})

test_that("compute_pca only returns available components", {
  dds <- suppressMessages(DESeq(makeExampleDESeqDataSet(n = 500, m = 4)))
  tr <- transform_counts(dds, "vst")
  p <- compute_pca(SummarizedExperiment::assay(tr), as.data.frame(colData(dds)))
  expect_equal(p$n_pc, 4)
  expect_length(p$percent, 4)
})

test_that("classify_genes handles NA padj", {
  st <- classify_genes(c(0.01, 0.01, 0.5, NA), c(2, -2, 3, 5), 0.05, 1)
  expect_equal(as.character(st), c("Up", "Down", "Not significant", "Not tested (padj = NA)"))
})

test_that("write_tsv writes an aligned header", {
  f <- tempfile()
  write_tsv(data.frame(gene = c("a", "b"), x = 1:2), f)
  expect_equal(readLines(f)[1], "gene\tx")
})
