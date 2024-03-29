context("Admixture analysis")

test_that("learn*partition behaves like on Matlab", {
  # Test data
  p1 <- c(0, .5, 1, 1.5)
  p2 <- c(seq(0, .5, .1), 1, 1, 1, 2)
  p3 <- c(.1, .1, .1, .5, .5, .5, 1, 1, 1)
  p4 <- c(.7, 1, 1, 1)

  # Testing learn_simple_partition
  expect_equal(
    object   = learn_simple_partition(p1, 2),
    expected = matrix(c(1, 1, 2, 2))
  )
  expect_equal(
    object   = learn_simple_partition(p2, 2),
    expected = matrix(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2))
  )
  expect_equal(
    object   = learn_simple_partition(p3, .5),
    expected = matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
  )
  expect_equal(
    object = learn_simple_partition(p4, 5),
    expected = matrix(c(1, 1, 1, 1))
  )
  expect_equal(
    object = learn_simple_partition(p4, .1),
    expected = matrix(c(1, 2, 2, 2))
  )

  # Testing learn_partition_modified
  expect_equal(
    object = learn_partition_modified(p4),
    expected = matrix(c(1, 2, 2, 2))
  )
})

test_that("type convertions behave like on Matlab", {
  expect_equal(ownNum2Str(1), "1")
  expect_equal(ownNum2Str(-123456789), "-123456789")
  expect_equal(ownNum2Str(0), "0")
  expect_error(ownNum2Str("a"))
  expect_equal(proportion2str(1), "1.00")
  expect_equal(proportion2str(0), "0.00")
  expect_equal(proportion2str(0.4), "0.40")
  expect_equal(proportion2str(0.89), "0.89")
  expect_equal(proportion2str(-0.4), "0.0-40") # also bugged in original
})

test_that("computeRows behaves like on Matlab", {
  # Matrices
  X <- matrix(1:9, 3, byrow = TRUE)
  Y <- matrix(9:1, 3, byrow = TRUE)
  Z <- matrix(c(-8, 2, -4, 0), byrow = TRUE)
  expect_equal(
    object   = computeRows(1, X, 3),
    expected = matrix(c(1, 4, 7))
  )
  expect_equal(
    object   = computeRows(2, X, 3),
    expected = matrix(c(1, 2, 7, 8, 13, 14))
  )
  expect_equal(
    object   = computeRows(10, X, 3),
    expected = matrix(c(1:10, 31:40, 61:70))
  )
  expect_equal(
    object   = computeRows(100, X, 3),
    expected = matrix(c(1:100, 301:400, 601:700))
  )
  expect_equal(
    object   = computeRows(1, Y, 3),
    expected = matrix(c(9, 6, 3))
  )
  expect_equal(
    object   = computeRows(2, Y, 3),
    expected = matrix(c(17, 18, 11, 12, 5, 6))
  )
  expect_equal(
    object   = computeRows(10, Y, 3),
    expected = matrix(c(81:90, 51:60, 21:30))
  )
  expect_equal(
    object   = computeRows(1, Z, 0),
    expected = matrix(, 1, 0)
  )
  expect_equal(
    object   = computeRows(1, Z, 5),
    expected = matrix(rep(-8, 5))
  )
  expect_equal(
    object   = computeRows(2, Z, 1),
    expected = matrix(rep(c(-17, -16), 1))
  )
  expect_equal(
    object   = computeRows(2, Z, 3),
    expected = matrix(rep(c(-17, -16), 3))
  )
  expect_equal(
    object   = computeRows(3, Z, 1),
    expected = matrix(rep(-26:-24, 1))
  )
  expect_equal(
    object   = computeRows(3, Z, 10),
    expected = matrix(rep(-26:-24, 10))
  )
})

test_that("computeIndLogml works like on Matlab", {
  expect_equivalent(computeIndLogml(10, 1), 2.3026, tol = .0001)
  expect_equivalent(computeIndLogml(0, 1), -Inf)
  expect_equivalent(computeIndLogml(1, 0), -Inf)
  expect_equivalent(computeIndLogml(0, 0), -Inf)
  expect_equivalent(computeIndLogml(-pi, -8), 3.2242, tol = .0001)
  expect_equivalent(computeIndLogml(2:3, 2), 2.3026, tol = .0001)
  expect_equivalent(computeIndLogml(matrix(8:5, 2), 100), 14.316, tol = .001)
  expect_equivalent(
    object    = computeIndLogml(matrix(8:5, 2), matrix(c(1, 3), 1)),
    expected  = 6.4118,
    tol       = .001
  )
  expect_equivalent(
    object    = computeIndLogml(matrix(8:5, 1), matrix(c(1, 3), 1)),
    expected  = 12.9717,
    tol       = .001
  )
  expect_equivalent(
    object = computeIndLogml(c(8, 1), c(-1.6, 5)),
    expected = complex(real = 6.4739, imaginary = pi),
    tol = .001
  )
})

test_that("suoritaMuutos works like on Matlab", {
  mx1 <- c(10, 5, 8)
  mx2 <- matrix(c(10, 9, 5, 8, 8, -7), 2)
  expect_equal(suoritaMuutos(10, 3, 1), 10)
  expect_equal(suoritaMuutos(mx1, 3, 1), c(10, 5, 8))
  expect_equal(suoritaMuutos(mx1, 3, 2), c(7, 8, 8))
  expect_equal(suoritaMuutos(mx1, 3, 3), c(7, 5, 11))
  expect_equal(suoritaMuutos(mx1, 2, 3), c(8, 5, 10))
  expect_equal(suoritaMuutos(mx1, -7, 3), c(17, 5, 1))
  expect_equal(suoritaMuutos(mx2, 0, 5), mx2)
  expect_equal(suoritaMuutos(mx2, 0, 5), mx2)
  expect_equal(suoritaMuutos(mx2, -3, 6), matrix(c(13, 9, 5, 8, 8, -10), 2))
})

test_that("laskeMuutokset4 works like on Matlab", {
  x <- admix1_muutokset$new()
  expect_equivalent(
    object = x$laskeMuutokset4(2, t(c(.4, 7)), c(8, 2), 3),
    expected = t(c(0, .3742)),
    tol = .0001
  )
})

test_that("etsiParas works like on Matlab", {
  mx1 <- t(c(.4, 7))
  expect_equal(etsiParas(2, mx1, c(8, 1), 8), c(.4, 7, 8))
  expect_equivalent(etsiParas(2, mx1, c(8, 1), 1), c(-1.6, 9, 3.1864), .0001)
  expect_equivalent(
    object = etsiParas(5, mx1, c(8, 1), -pi),
    expected = c(-4.6, 12, 3.8111),
    tol = .001
  )
  expect_equivalent(
    object = etsiParas(-.5, mx1, c(-1, 0), -10),
    expected = c(7.4, 0, complex(real = 1.8563, imaginary = 3.1416)),
    tol = .0001
  )
})

test_that("computePersonalAllFreqs works like on Matlab", {
  expect_equal(computePersonalAllFreqs(1, 1:4, c(15, 5, 10, 40), 1), 15)
  mx <- matrix(c(15, 10, 5, 40), 2)
  expect_equal(computePersonalAllFreqs(1, 1:4, mx, 1), c(15, 40))
  expect_equal(computePersonalAllFreqs(1, 1:3, mx, 1), c(15, 40))
  expect_equal(computePersonalAllFreqs(1, 1:2, mx, 1), c(15, 40))
})

test_that("simuloiAlleeli works like on Matlab", {
  sk <- 2
  vk <- 1:3
  ra <- array(1:12, c(2, 2, 3))
  mx1 <- matrix(c(3, 5, 0, 9), 2)
  mx2 <- matrix(c(3, 5, 0, 9, 5, 8), 2)
  expect_equal(simuloiAlleeli(sk, 1, 1), 1)
  expect_equal(simuloiAlleeli(vk, 1, 2), 1)
  expect_equal(simuloiAlleeli(ra, 2, 1), 1)
  expect_equal(simuloiAlleeli(mx1, 1, 2), 2)
  expect_equal(simuloiAlleeli(mx2, 1, 3), 1)
})

test_that("simulateIndividuals works like on Matlab", {
  set.seed(2)
  expect_equal(
    object = simulateIndividuals(1, 3, 2, 0, .2),
    expected = matrix(c(1, -999, 1), ncol = 1)
  )
  expect_equal(
    object = simulateIndividuals(5, 3, 1:3, 4, 0),
    expected = matrix(rep(-999, 15 * 3), 15)
  )
  expect_equal(
    object = simulateIndividuals(3, 3, 2, 1, 1),
    expected = matrix(rep(1, 9), 9)
  )
  set.seed(2)
  expect_equal(
    object = sum(simulateIndividuals(3, 3, 2, 1, .5) == 1),
    expected = 6
  )
})

test_that("simulateAllFreqs works as expected", {
  empty_mt <- matrix(NA, 0, 0)
  expect_equivalent(suppressWarnings(simulateAllFreqs(3)), empty_mt)
  expect_equivalent(suppressWarnings(simulateAllFreqs(3:5)), empty_mt)
  expect_equivalent(
    object = suppressWarnings(simulateAllFreqs(matrix(1:4, 2))),
    expected = empty_mt
  )
})

test_that("computeAllFreqs2 works as expected", {
  expect_equivalent(computeAllFreqs2(10), matrix(NA, 0, 0))
})

test_that("poistaLiianPienet works as expected", {
  expect_equal(poistaLiianPienet(100, matrix(1:4, 2), 0), 100)
  expect_equal(poistaLiianPienet(100, matrix(1:4, 2), -5), 100)
})

test_that("noIndex works properly", {
  abcd_vec <- letters[1:4]
  abcd_mat <- matrix(abcd_vec, 2)
  abcdef_mat <- matrix(letters[1:6], 2)
  efg_vec <- letters[5:7]
  expect_equal(noIndex(abcd_vec, 1:6), abcd_vec)
  expect_equal(noIndex(abcd_vec, 1:3), abcd_vec[-4])
  expect_equal(noIndex(abcd_vec, 1:2), abcd_vec)
  expect_equal(noIndex(abcd_vec, efg_vec), abcd_vec[-4])
  expect_equal(noIndex(abcd_mat, 1), abcd_mat[, 1])
  expect_equal(noIndex(abcd_mat, 2), abcd_mat[, 1])
  expect_equal(noIndex(abcdef_mat, 1:2), abcdef_mat[, 1:2])
  expect_equal(noIndex(abcdef_mat, abcd_mat), abcdef_mat[, 1:2])
})
