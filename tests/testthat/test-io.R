test_that("read_counts returns an integer matrix", {
  f <- write_tmp(c("gene\ts1\ts2", "g1\t1\t2", "g2\t3\t4"))
  x <- read_counts(f)
  expect_true(is.integer(x$counts))
  expect_equal(dim(x$counts), c(2L, 2L))
  expect_length(x$messages, 0)
})

test_that("read_counts rounds non-integer counts and warns", {
  f <- write_tmp(c("gene\ts1\ts2", "g1\t1.4\t2.6", "g2\t3\t4"))
  x <- read_counts(f)
  expect_equal(unname(x$counts[1, ]), c(1L, 3L))
  expect_match(x$messages, "rounded")
})

test_that("read_counts rejects invalid matrices", {
  expect_error(read_counts(write_tmp(c("gene\ts1\ts2", "g1\t1\tNA", "g2\t3\t4"))), "missing")
  expect_error(read_counts(write_tmp(c("gene\ts1\ts2", "g1\t1\t-2"))), "negative")
  expect_error(read_counts(write_tmp(c("gene\ts1\ts2", "g1\t1\tabc"))), "Non-numeric")
  expect_error(read_counts(write_tmp(c("gene,s1,s2", "g1,1,2"))), "at least two")
})

test_that("read_design uses the first column as IDs and any column name", {
  f <- write_tmp(c("id\tcondition\tbatch", "s1\tctrl\tA", "s2\ttrt\tB"))
  x <- read_design(f)
  expect_equal(rownames(x$design), c("s1", "s2"))
  expect_true(all(vapply(x$design, is.factor, logical(1))))
  # level order follows first appearance, so the first level is the reference
  expect_equal(levels(x$design$condition), c("ctrl", "trt"))
})

test_that("read_design sanitises column names and rejects duplicated IDs", {
  x <- read_design(write_tmp(c("id\tcell type", "s1\tA", "s2\tB")))
  expect_equal(names(x$design), "cell.type")
  expect_match(x$messages, "syntactically valid")
  expect_error(read_design(write_tmp(c("id\tc", "s1\tA", "s1\tB"))), "Duplicated")
})

test_that("match_samples reorders and reports mismatches both ways", {
  counts <- matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s2", "s1")))
  design <- data.frame(c = factor(c("a", "b")), row.names = c("s1", "s2"))
  expect_equal(colnames(match_samples(counts, design)), c("s1", "s2"))

  design_extra <- data.frame(c = factor(c("a", "b", "b")), row.names = c("s1", "s2", "s3"))
  expect_error(match_samples(counts, design_extra), "Only in design: s3")
  design_missing <- data.frame(c = factor("a"), row.names = "s1")
  expect_error(match_samples(counts, design_missing), "Only in counts: s2")
})

test_that("build_formula puts the variable of interest last", {
  expect_equal(deparse(build_formula("condition")), "~condition")
  expect_equal(deparse(build_formula("condition", c("batch", "condition"))), "~batch + condition")
})
