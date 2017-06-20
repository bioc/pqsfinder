context("Regressions")

expect_equal_pv_vectors <- function(pv_a, pv_b) {
  expect_equal(maxScores(pv_a), maxScores(pv_b))
  expect_equal(density(pv_a), density(pv_b))
}

expect_equal_pv_coords <- function(pv_a, pv_b) {
  expect_equal(length(pv_a), length(pv_b))
  expect_equal(start(pv_a), start(pv_b))
  expect_equal(end(pv_a), end(pv_b))
  expect_equal(score(pv_a), score(pv_b))
}

expect_no_overlaps <- function(pv) {
  cnts <- numeric(length(subject(pv)))
  st <- start(pv)
  ed <- end(pv)
  for (i in 1:length(pv)) {
    cnts[st[i]:ed[i]] <- cnts[st[i]:ed[i]] + 1
  }
  expect_equal(sum(cnts > 1), 0)
}

test_that("the result on test seq is the same as gives pqsfinder-1.4.4-patched", {
  load("pqsfinder_1_4_4_patched.RData")
  
  pv_d <- pqsfinder(test_seq, strand = "+")
  pv_r <- pqsfinder(test_seq, strand = "+", run_re = "G{1,10}.{0,10}G{1,10}")
  
  print(test_seq)
  library(rtracklayer)
  export(as(pv_d, "GRanges"), "pv_d.gff", version = "3")
  cat(readLines("pv_d.gff"), sep = "\n")
  
  cat("pv_d, pv_r\n")
  expect_equal_pv_vectors(pv_d, pv_r)
  expect_no_overlaps(pv_d)
  expect_no_overlaps(pv_r)
  
  cat("pqsfinder_1_4_4_d, pqsfinder_1_4_4_r\n")
  expect_equal_pv_vectors(pqsfinder_1_4_4_patched_d, pqsfinder_1_4_4_patched_r)
  expect_no_overlaps(pqsfinder_1_4_4_patched_d)
  expect_no_overlaps(pqsfinder_1_4_4_patched_r)
  
  cat("pv_d, pqsfinder_1_4_4_d\n")
  pv_i <- pv_d[start(pv_d) %in% start(pqsfinder_1_4_4_patched_d)]
  expect_equal_pv_coords(pv_i, pqsfinder_1_4_4_patched_d)
})
