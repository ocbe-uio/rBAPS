context("Basic Matlab functions")

test_that("rand works properly", {
    expect_equal(dim(rand()), c(1, 1))
    expect_equal(dim(rand(1, 2)), c(1, 2))
    expect_equal(dim(rand(3, 2)), c(3, 2))
})

test_that("repmat works properly", {
    mx <- matrix(1:4)
    expect_equal(repmat(1:4, 1), mx)
    expect_equal(repmat(1:4, 2), t(cbind(rbind(mx, mx), rbind(mx, mx))))
    expect_equal(
        object = repmat(1:4, c(2, 3)),
        expected = t(cbind(rbind(mx, mx), rbind(mx, mx), rbind(mx, mx)))
    )
    expect_equal(
        object = repmat(1:4, c(4, 1)),
        expected = rbind(mx, mx, mx, mx)
    )
})

test_that("zeros and ones work as expected", {
    expect_equal(zeros(1), matrix(0, 1))
    expect_equal(zeros(2), matrix(0, 2, 2))
    expect_equal(zeros(2, 1), matrix(0, 2, 1))
    expect_equal(zeros(1, 10), matrix(0, 1, 10))
    expect_equal(ones(8), matrix(1, 8, 8))
    expect_equal(ones(5, 2), matrix(1, 5, 2))
    expect_equal(ones(2, 100), matrix(1, 2, 100))
})

test_that("times works as expected", {
    expect_equal(times(9, 6), as.matrix(54))
    expect_equal(times(9, c(.8, 9)), as.matrix(c(7.2, 81)))
    expect_equal(times(c(.8, 9), 5), as.matrix(c(4, 45)))
    expect_equal(times(matrix(1:4, 2), 5), matrix(c(5, 10, 15, 20), 2))
    expect_equal(times(5, matrix(1:4, 2)), matrix(c(5, 10, 15, 20), 2))
    expect_equal(times(matrix(1:4, 2), c(10, 3)), matrix(c(10, 20, 9, 12), 2))
    expect_equal(times(c(10, 3), matrix(1:4, 2)), matrix(c(10, 20, 9, 12), 2))
    expect_equal(
        object = times(matrix(c(10, -5, 3, 9), 2), matrix(1:4, 2)),
        expected = matrix(c(10, -10, 9, 36), 2)
    )
})