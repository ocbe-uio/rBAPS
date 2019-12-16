context("Admixture analysis")


test_that("learn_simple_partition behaves like Matlab", {
    p1 <- c(0, .5, 1, 1.5)
    p2 <- c(seq(0, .5, .1), 1, 1, 1, 2)
    p3 <- c(.1, .1, .1, .5, .5, .5, 1, 1, 1)
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
})
