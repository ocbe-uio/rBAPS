context("greedyPopMix functions")

test_that("Auxiliary functions work properly", {
  x <- matrix(11:16, 3)
  y <- matrix(2:7, 3)
  z <- list(
    popnames2 = matrix(c(11:13, seq(1.5, 2.5, 0.5)), 3),
    rowsFromInd = 2
  )
  x2 <- matrix(seq(4, 14, 2), 3)
  expect_equal(findOutRowsFromInd(x, y, "Diploid"), z)
  expect_equal(
    getPopDistancesByKL(x2),
    list(
      Z = matrix(c(c(1, 101:198), c(2:100), rep(0, 99)), nrow = 99, ncol = 3),
      distances = as.matrix(rep(0, 4950))
    )
  )
  x3 <- matrix(c(9, 4, 1, 8, 5, 2, 1, 2, 3), 3)
  y3 <- list(
    newData = matrix(c(3:1, 3:1, 1:3), 3),
    rows = matrix(c(1:3, 1:3), 3),
    alleleCodes = matrix(c(1, 4, 9, 2, 5, 8), 3),
    noalle = c(3, 3),
    adjprior = matrix(rep(3/9, 6), 3),
    priorTerm = 5.9125
  )
  expect_equal(handlePopData(x3), y3)
})
