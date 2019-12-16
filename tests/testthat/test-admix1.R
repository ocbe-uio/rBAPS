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

test_that("ownNum2Str behaves like on Matlab", {
    expect_equal(ownNum2Str(1), "1")
    expect_equal(ownNum2Str(-123456789), "-123456789")
    expect_equal(ownNum2Str(0), "0")
    expect_error(ownNum2Str("a"))
})