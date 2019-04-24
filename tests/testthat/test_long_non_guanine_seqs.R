context("Long non-guanine sequences")

test_that("nothing found in long non-guanine sequences (should be fast)", {
  expect_equal(length(pqsfinder(DNAString(strrep("TA", times = 2000000)), min_score = 26)), 0)
})
