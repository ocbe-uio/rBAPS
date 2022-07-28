context("Spatial mixture")

test_that("functions work with basic input", {
  q <- c(1, 3, 0, 2)
  w <- array(c(1, 3, 2, 4, 5, 7, 6, 8), c(2, 2, 2))
  e <- computeCounts(q, 4, 5, w)
  r <- matrix(c(5, 3, 3, 6, 4, 4, 7, 9, 1), 3)
  expect_equal(e$cliqcounts, t(c(1, 1, 1, 0, 0)))
  expect_equal(e$sepcounts, t(c(0, 0, 0, 1, 0)))
  expect_equal(computeDiffInCliqCounts(r, 4), matrix(c(0, 1, 1)))
  expect_equal(computeDiffInCliqCounts(r, 3), matrix(c(0, 1, 1)))
  expect_equal(computeDiffInCliqCounts(r, 5), matrix(c(1, 0, 0)))
  expect_equal(computeDiffInCliqCounts(r, 0), matrix(c(0, 0, 0)))
})
