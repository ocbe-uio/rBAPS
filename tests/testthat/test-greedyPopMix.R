context("greedyPopMix functions")

test_that("Auxiliary functions work properly", {
  # Testing findOutRowsFromInd -------------------------------------------------
  x <- matrix(11:16, 3)
  y <- matrix(2:7, 3)
  z <- list(
    popnames2 = matrix(c(11:13, seq(1.5, 2.5, 0.5)), 3),
    rowsFromInd = 2
  )
  expect_equal(findOutRowsFromInd(x, y, "Diploid"), z)

  # Testing getPopDistancesByKL ------------------------------------------------
  set.seed(5562143)
  globals$COUNTS <- array(rpois(10 ^ 3, lambda = 10), dim = c(10, 10, 10))
  x2 <- matrix(seq(4, 14, 2), 3)
  Z <- matrix(
    c(
      5, 8, 0.004410383, 11, 10, 0.009386554, 2, 7, 0.010046641,
      3, 12, 0.011518562, 1, 4, 0.013682110, 14, 6, 0.027109229,
      15, 13, 0.032833011, 16, 9, 0.045394132, 17, 18, 0.086729590
    ),
    nrow = 9, ncol = 3, byrow = TRUE
  )
  distances <- matrix(
    c(
      0.015767099, 0.020110943, 0.013682110, 0.010214691, 0.050375802,
      0.019692440, 0.010294099, 0.053740306, 0.028410233, 0.034496635,
      0.025356695, 0.015267764, 0.045409885, 0.010046641, 0.012797950,
      0.065436956, 0.026308651, 0.019896390, 0.006193117, 0.011735103,
      0.036077638, 0.011518562, 0.029233839, 0.010798344, 0.013262098,
      0.034671122, 0.032833011, 0.024732472, 0.033976345, 0.040243817,
      0.020374790, 0.024770297, 0.004410383, 0.024329842, 0.009386554,
      0.052217370, 0.027109229, 0.038127778, 0.016888110, 0.024061302,
      0.086729590, 0.038787638, 0.045394132, 0.005206076, 0.042606915
    ),
    nrow = 45, ncol = 1,
    byrow = TRUE
  )
  expect_equal(getPopDistancesByKL(x2), list(Z = Z, distances = distances))

  # Testing handlePopData ------------------------------------------------------
  x3 <- matrix(c(9, 4, 1, 8, 5, 2, 1, 2, 3), 3)
  y3 <- list(
    newData     = matrix(c(3: 1, 3: 1, 1: 3), 3),
    rows        = matrix(c(1: 3, 1: 3), 3),
    alleleCodes = matrix(c(1, 4, 9, 2, 5, 8), 3),
    noalle      = matrix(c(3, 3), 1),
    adjprior    = matrix(rep(3 / 9, 6), 3),
    priorTerm   = 5.9125
  )
  expect_equal(handlePopData(x3), y3, tol = 1e-4)
})
