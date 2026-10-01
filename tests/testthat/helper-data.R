write_tmp <- function(lines) {
  f <- tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

example_counts <- function() file.path("..", "..", "inst", "extdata", "pasilla_counts.tsv")
example_design <- function() file.path("..", "..", "inst", "extdata", "pasilla_design.tsv")
