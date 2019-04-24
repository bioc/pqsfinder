context("Zero-length loop")

test_that("only one zero-length loop is allowed", {
  expect_equal(length(pqsfinder(DNAString("GGTGGGGGG"), min_score = 26)), 0)
  expect_equal(length(pqsfinder(DNAString("GGGGTGGGG"), min_score = 26)), 0)
  expect_equal(length(pqsfinder(DNAString("GGGGGGTGG"), min_score = 26)), 0)
  expect_equal(length(pqsfinder(DNAString("GGTGGTGGGG"), min_score = 26)), 1)
  expect_equal(length(pqsfinder(DNAString("GGTGGGGTGG"), min_score = 26)), 1)
  expect_equal(length(pqsfinder(DNAString("GGGGTGGTGG"), min_score = 26)), 1)
})
