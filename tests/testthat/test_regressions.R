context("Regressions")

expect_equal_pv <- function(pv_a, pv_b) {
  expect_equal(length(pv_a), length(pv_b))
  expect_equal(maxScores(pv_a), maxScores(pv_b))
  expect_equal(density(pv_a), density(pv_b))
  expect_equal(start(pv_a), start(pv_b))
  expect_equal(end(pv_a), end(pv_b))
}

test_that("the result on test seq is the same as gives pqsfinder-1.4.4", {
  load("pqsfinder_1_4_4.RData")
  
  pv_d <- pqsfinder(test_seq, strand = "+")
  pv_r <- pqsfinder(test_seq, strand = "+", run_re = "G{1,10}.{0,10}G{1,10}")
  
  expect_equal_pv(pv_d, pv_r)
  expect_equal_pv(pqsfinder_1_4_4_d, pqsfinder_1_4_4_r)
  
  expect_equal_pv(pv_d, pqsfinder_1_4_4_d)
  expect_equal_pv(pv_r, pqsfinder_1_4_4_r)
})
