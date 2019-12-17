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
    expect_equal(proportion2str(-0.4), "0.0-40")  # also bugged in original
    # TODO: fix after release, as long as it doesn't break anything else
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