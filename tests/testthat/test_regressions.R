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

test_that("sequences pqs parts can be extracted", {
  
  library(stringr)
  
  test_seq <- DNAString("GGGTAGTGGTTTTGGGTTTGGGAAAAAAAAAAAAAAGGGTTTGGAGGAAATTTGGGGAGGGG")
  pv <- pqsfinder(test_seq, strand = "+")
  pv_m <- elementMetadata(pv)
  
  r1_s <- start(pv)
  r1_e <- r1_s + pv_m$rl1
  
  l1_s <- r1_e
  l1_e <- l1_s + pv_m$ll1
  
  r2_s <- l1_e
  r2_e <- r2_s + pv_m$rl2
  
  l2_s <- r2_e
  l2_e <- l2_s + pv_m$ll2
  
  r3_s <- l2_e
  r3_e <- r3_s + pv_m$rl3
  
  l3_s <- r3_e
  l3_e <- l3_s + pv_m$ll3
  
  r4_s <- l3_e
  r4_e <- end(pv)
  
  # apply end point corrections
  r1_e <- r1_e - 1
  l1_e <- l1_e - 1
  r2_e <- r2_e - 1
  l2_e <- l2_e - 1
  r3_e <- r3_e - 1
  l3_e <- l3_e - 1
  
  r1 <- str_sub(test_seq, r1_s, r1_e)
  l1 <- str_sub(test_seq, l1_s, l1_e)
  r2 <- str_sub(test_seq, r2_s, r2_e)
  l2 <- str_sub(test_seq, l2_s, l2_e)
  r3 <- str_sub(test_seq, r3_s, r3_e)
  l3 <- str_sub(test_seq, l3_s, l3_e)
  r4 <- str_sub(test_seq, r4_s, r4_e)
  
  pqs_seqs = sprintf("[%s]%s[%s]%s[%s]%s[%s]", r1, l1, r2, l2, r3, l3, r4)
  
  expect_equal(nchar(pqs_seqs), width(pv) + 8)
})
