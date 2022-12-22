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
  # TODO: test elim_order()
  # TODO: test triangulate()
  # TODO: test myintersect()
  # TODO: test findCliques()
  # TODO: test cliques_to_jtree()
  # TODO: test lakseKlitik()
})
