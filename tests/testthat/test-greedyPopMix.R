context("greedyPopMix functions")

test_that("Auxiliary functions work properly", {
  x <- matrix(11:16, 3)
  y <- matrix(2:7, 3)
  z <- list(
    popnames2 = matrix(c(11:13, seq(1.5, 2.5, 0.5)), 3),
    rowsFromInd = 2
  )
  expect_equal(findOutRowsFromInd(x, y, "Diploid"), z)
})
