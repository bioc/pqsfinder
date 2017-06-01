context("Zero-length loop")

test_that("only one zero-length loop is allowed", {
  expect_equal(length(pqsfinder(DNAString("GGTGGGGGG"))), 0)
  expect_equal(length(pqsfinder(DNAString("GGGGTGGGG"))), 0)
  expect_equal(length(pqsfinder(DNAString("GGGGGGTGG"))), 0)
  expect_equal(length(pqsfinder(DNAString("GGTGGTGGGG"))), 1)
  expect_equal(length(pqsfinder(DNAString("GGTGGGGTGG"))), 1)
  expect_equal(length(pqsfinder(DNAString("GGGGTGGTGG"))), 1)
})
