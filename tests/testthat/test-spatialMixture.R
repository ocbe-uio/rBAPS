context("Spatial mixture")

test_that("functions work with basic input", {
  q <- c(1, 3, 0, 2)
  w <- array(c(1, 3, 2, 4, 5, 7, 6, 8), c(2, 2, 2))
  e <- computeCounts(q, 4, 5, w)
  r <- matrix(c(5, 3, 3, 6, 4, 4, 7, 9, 1), 3)
  t <- 1:5
  y <- c(4, 3, 1, 4, 8)
  expect_equal(e$cliqcounts, t(c(1, 1, 1, 0, 0)))
  expect_equal(e$sepcounts, t(c(0, 0, 0, 1, 0)))
  expect_equal(computeDiffInCliqCounts(r, 4), matrix(c(0, 1, 1)))
  expect_equal(computeDiffInCliqCounts(r, 3), matrix(c(0, 1, 1)))
  expect_equal(computeDiffInCliqCounts(r, 5), matrix(c(1, 0, 0)))
  expect_equal(computeDiffInCliqCounts(r, 0), matrix(c(0, 0, 0)))
  expect_equal(mysetdiff(t, y), c(2, 5))
})

test_that("testaaKoordinaatit works as expected", {
  m1 <- matrix(c(11.1, 22.2, 33.3, 44.4), 2, byrow = TRUE)
  m2 <- matrix(c(11.1, 22.2, 33.3, 44.4, 55.5, 66.6), byrow = TRUE)
  m3 <- matrix(c(11, 22.2, 11, 44.4, 11, 66.6), ncol = 2, byrow = TRUE)
  expect_equal(
    testaaKoordinaatit(2, m1), list("viallinen" = 0, "coordinates" = m1)
  )
  expect_warning(testaaKoordinaatit(2, m2), "Wrong coordinates dimension!")
  expect_equal(
    {set.seed(5676402); testaaKoordinaatit(3, m3, FALSE)},
    list(
      "viallinen" = 0,
      "coordinates" = matrix(
        c(11.99, 22.2, 11.50, 44.4, 11.00, 66.6), ncol = 2, byrow = TRUE)
      )
  )
})

test_that("lakseKlitik() and subfunctions produce expected output", {
  expect_equal(neighbors(matrix(c(11, 22, 33, 44), 2), 2), c(1, 2))
  expect_equal(myintersect(matrix(1:4, 2), matrix(2:5, 2)), 2:4)
  expect_equal(myintersect(matrix(1:4, 2), matrix(5:8, 2)), integer(0))
  expect_equal(myintersect(matrix(1:4, 2), matrix(4:7, 2)), 4)
  expect_true(myisvector(runif(1)))
  expect_true(myisvector(matrix(runif(1))))
  expect_true(myisvector(runif(2)))
  expect_true(myisvector(matrix(runif(2))))
  expect_true(myisvector(rand(2, 1)))
  expect_true(myisvector(rand(1, 2)))
  expect_false(myisvector(rand(2, 2)))
  expect_equal(mysize(rand(1, 1)), 1)
  expect_equal(mysize(rand(2, 1)), 2)
  expect_equal(mysize(rand(1, 2)), 2)
  expect_equal(mysize(rand(2, 2)), c(2, 2))
  expect_equal(dec2bitv(1, 2), c(0, 1))
  expect_equal(dec2bitv(5, 2), c(1, 0, 1))
  expect_equal(dec2bitv(5, 5), c(0, 0, 1, 0, 1))
  expect_equal(dec2bitv(5, 1), c(1, 0, 1))
  expect_equal(dec2bitv(5, 0), c(1, 0, 1))
  expect_equal(dec2bitv(10, 1), c(1, 0, 1, 0))
  expect_equal(dec2bitv(10, 5), c(0, 1, 0, 1, 0))
  expect_equal(dec2bitv(10, 10), c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0))
  # TODO: test ind2subv()
  # TODO: test argmin()
  # TODO: test elim_order()
  # TODO: test triangulate()
  # TODO: test mysubset()
  # TODO: test findCliques()
  # TODO: test cliques_to_jtree()
  # TODO: test minimum_spanning_tree()
  # TODO: test myunion()
  # TODO: ... and anythin left from findCliques.m
  # TODO: test lakseKlitik()
})

test_that("testFastaData() produces same output as on MATLAB", {
  msa <- system.file("ext", "seqs.fa", package = "rBAPS")
  test_msa <- testFastaData(msa)
  expect_equal(test_msa$ninds, 515)
  expect_equal(dim(test_msa$data), c(515, 745))
  expect_named(table(test_msa$data), c("-9", as.character(1:515)))
})
