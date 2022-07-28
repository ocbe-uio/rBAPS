context("Spatial mixture")

test_that("functions work with basic input", {
  x <- c(1, 3, 0, 2)
  y <- array(c(1, 3, 2, 4, 5, 7, 6, 8), c(2, 2, 2))
  z <- computeCounts(x, 4, 5, y)
  expect_equal(z$cliqcounts, t(c(1, 1, 1, 0, 0)))
  expect_equal(z$sepcounts, t(c(0, 0, 0, 1, 0)))
})
