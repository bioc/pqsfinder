context("Exhaustive search")

test_that("expansion of sequence won't prevent to find at least the same solution", {
  a <- DNAString("CCCCCCGGGTGGGTGAGGTGGGAA")
  d <- DNAString("CCCCCCGGGTGGGTGAGGTGGGAAAAGGGCATTATTCTAGGGGAAAAAAAGGGAAAGGGAAA")
  
  pqs_a <- pqsfinder(a, min_score = 26)
  expect_equal(start(pqs_a), 7)
  expect_equal(width(pqs_a), 16)
  
  pqs_d <- pqsfinder(d, min_score = 26)
  expect_equal(start(pqs_d), 7)
  expect_true(width(pqs_d) >= width(pqs_a))
  expect_true(score(pqs_d) >= score(pqs_a))
})
