context("Regressions")
library(stringi)

pqs_max_scores <- function(pv, i) {
  maxScores(pv)[start(pv)[i]:end(pv)[i]]
}

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

# only works for single strand search outputs
expect_no_overlaps <- function(pv) {
  cnts <- numeric(length(subject(pv)))
  st <- start(pv)
  ed <- end(pv)
  for (i in 1:length(pv)) {
    cnts[st[i]:ed[i]] <- cnts[st[i]:ed[i]] + 1
  }
  expect_equal(sum(cnts > 1), 0)
}

test_that("pqsfinder-1.4.4-patched results are consistent", {
  load("pqsfinder_1_4_4_patched.RData")
  
  cat(" compare pqsfinder_1_4_4_d, pqsfinder_1_4_4_r\n")
  expect_equal_pv_vectors(pqsfinder_1_4_4_patched_d, pqsfinder_1_4_4_patched_r)
  expect_no_overlaps(pqsfinder_1_4_4_patched_d)
  expect_no_overlaps(pqsfinder_1_4_4_patched_r)
})

test_that("the result on test seq is the same as gives pqsfinder-1.4.4-patched", {
  load("pqsfinder_1_4_4_patched.RData")
  
  cat(" run default pqsfinder\n")
  pv_d <- pqsfinder(test_seq, strand = "+", deep = FALSE)
  cat(" run pqsfinder using boost regex engine\n")
  pv_r <- pqsfinder(test_seq, strand = "+", run_re = "G{1,10}.{0,10}G{1,10}", deep = FALSE)
  
  cat(" compare pv_d, pv_r\n")
  expect_no_overlaps(pv_d)
  expect_no_overlaps(pv_r)
  
  # reflect slight change in the mechanism how same-scoring PQS overlaps are resolved
  start(pqsfinder_1_4_4_patched_d)[14] <- 14962
  width(pqsfinder_1_4_4_patched_d)[14] <- 29
  
  cat(" compare pv_d, pqsfinder_1_4_4_d\n")
  pv_i <- pv_d[start(pv_d) %in% start(pqsfinder_1_4_4_patched_d)]
  expect_equal_pv_coords(pv_i, pqsfinder_1_4_4_patched_d)
})

test_that("fast computation give same pqs as original slow one on real test_seq", {
  load("pqsfinder_1_4_4_patched.RData")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  pv_fast_p <- pv_fast[strand(pv_fast) == "+"]
  pv_fast_m <- pv_fast[strand(pv_fast) == "-"]
  
  pv_slow_p <- pv_slow[strand(pv_slow) == "+"]
  pv_slow_m <- pv_slow[strand(pv_slow) == "-"]
  
  expect_no_overlaps(pv_fast_p)
  expect_no_overlaps(pv_fast_m)
  expect_no_overlaps(pv_slow_p)
  expect_no_overlaps(pv_slow_m)
  
  # expect_equal_pv_vectors(pv_slow, pqsfinder_1_4_4_patched_d)
  
  pv_i <- pv_fast[start(pv_fast) %in% start(pv_slow)]
  expect_equal_pv_coords(pv_i, pv_slow)
})

test_that("there no overshadowing", {
  test_seq <- DNAString("GGGCATGGCCCACCAGGAGGGTGGCTCGGGTGGGGGACAAGCTGGTAGGCAGGGCCAAGGAGCTGAGAGGCTACACGGGAGGGAGCTGACCCACAGGACCAGGACAGGGGGCTTGGAGGAGGCGAGCAGAGGAGCTGGGG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("maxScores vector is a max of maxScores on sense and antisense", {
  test_seq <- DNAString(stri_rand_shuffle(strrep("ACCGGT", 500)))
  pv <- pqsfinder(test_seq, deep = TRUE)
  pv_p <- pqsfinder(test_seq, strand = "+", deep = TRUE)
  pv_m <- pqsfinder(test_seq, strand = "-", deep = TRUE)
  expect_equal(maxScores(pv), pmax(maxScores(pv_p), maxScores(pv_m)))
})

test_that("density vector is a sum of density on sense and antisense", {
  test_seq <- DNAString(stri_rand_shuffle(strrep("ACCGGT", 500)))
  pv <- pqsfinder(test_seq, deep = TRUE)
  pv_p <- pqsfinder(test_seq, strand = "+", deep = TRUE)
  pv_m <- pqsfinder(test_seq, strand = "-", deep = TRUE)
  expect_equal(density(pv), density(pv_p) + density(pv_m))
})

test_that("results are the same for slow and fast computation on random string", {
  test_seq <- DNAString(stri_rand_shuffle(strrep("ACGGGT", 1000)))
  
  system.time(pv_fast <- pqsfinder(test_seq, deep = FALSE))
  system.time(pv_slow <- pqsfinder(test_seq, deep = TRUE))
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq1", {
  test_seq <- DNAString("GGGGGGGCAGCCGGCTCGTTTCGGAGAGCGGGCCGGGCAAGGGTGAACACACAGCGGGC")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq2", {
  test_seq <- DNAString("GGGGTGCCCTAGGGGAGGGGGTGGGGAGGGACCGCCGGGCGACGGAGGGGTCGGTGCTTAGGACTGTCGGGGCGAAACGGCAGCGG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq3", {
  test_seq <- DNAString("GTACTGGGCCGGGGGCGGGCGGTGGCTGGCTAGCTGGGGGGGGGAGGCGCCGGTCAACGGGGGGCCAAGGAGGAGGG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq4", {
  test_seq <- DNAString("GGCTAGCGGGTGGAGGGCTCGAGATGCGGAGTGTTGGTAGGTGGGAGAGCCTTAGATACGGGTAACGAGGGGTCGCGCGGCAGTGTTGAGACTTGTGGGCAGCGTGGTAGGGGCGCGGGCGAAGGGGGGCTGGGCGGGTGTCGGGGAAGGAGGGGGCG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq5", {
  test_seq <- DNAString("ATGAGTGGGGCGTTTGCGGGTGCGGGAAGATGGGAGGCGCATTCGGGACAGGAGGGGGACATTTCGGAACGTGAGGGCCGGGGGTGTAAGGGGACAGGGGGTGAGGGCGGGGAGGGGCTTGGTTCGAGTATCGTTGGGGGGTCGGCTGGGGCAGGAGCGATGGGGCGAACTCGTGGGCGGGCTTAGGTGCGAGGTGGATACGTCGAAGCTAGACGTAGGGGGACAGCGCGGGGGGGAGCGATAATGTTCTG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq6", {
  test_seq <- DNAString("GGATTGGGGGACCGGGGGTTGGTGGGGTGGATACTGGGGGGGGGGGGCAGTTCTCAGTTAAAGTGCGGTGCCGGAGGGCGACGGGGTGAATATGTAGACAGGGGGCGCCGGGGGGGGCATTGTGGGGCAGAATAGGGATGGGTATGATGTCGGGCG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq7", {
  test_seq <- DNAString("GGGGCGAGGCAGAGGTGGGCAAGGACAAGGGGCAGGCAGGAATTAAGAGTGAGAATGCGGCTGCAGCCAGGAAGTGAGGTCTGCGCGCAGCCAGGGACTGGAGCACTCCATCAGGCGAGCTGGAGAGAAGGGGGCTTGGAACCAGGATTCGGGAGGAGGAGCTGCTGCCAGGGCTGACAGGGTCAAGGGCATGGCCCACCAGGAGGGTGGCTCGGGTGGGGGACAAGCTGGTAGGCAG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq8", {
  test_seq <- DNAString("GGGCCAAGGAGCTGAGAGGCTACACGGGAGGGAGCTGACCCACAGGACCAGGACAGGGGGCTTGGAGGAGGCGAGCAGAGGAGCTGGGGCTTTGCACAGGGTGAGGAGAGAGGGCTCGACCATATGGGAAAGGGAGGAACCACCCATGGG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("results are correct on buggy seq9", {
  test_seq <- DNAString("GGGTACAGGGTGTGCTAGGGGGGCGGGGGTGCTCGGGGAGTTGGGTTGAGAAGTGAGTGGGGGATGGGTCACTGGCGGGGCGGGGTAGGGGGGG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
})

test_that("the shortest possible between same-scoring PQS is found", {
  test_seq <- DNAString("GGGGCCTGTCAGGGGGTCGGGGAAGGGGGGGCTAGGGGAG")
  
  pv_fast <- pqsfinder(test_seq, deep = FALSE)
  pv_slow <- pqsfinder(test_seq, deep = TRUE)
  
  expect_equal(width(pv_fast), 29)
  expect_equal(width(pv_slow), 29)
  
  expect_equal_pv_coords(pv_fast, pv_slow)
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
