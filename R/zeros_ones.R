#' @title Matrix of zeros or ones
#' @description Generates a square or rectangular matrix of zeros or ones
#' @param n scalar or 2D vector
#' @param x value to fill matrix with
#' @return n-by-n matrix filled with `x`
#' @details This is a wrapper function to replicate the behavior of the
#' `zeros()` and the `ones()` functions on Matlab
#' @note Actually works for any `x`, but there's no need to bother imposing
#' validation controls here.
zeros_or_ones <- function(n, x) {
    # Expanding n to length 2 if necessary
    if (length(n) == 1) n <- c(n, n)

    # Returning a matrix or an array
    if (length(n) == 2) {
        return(matrix(x, n[1], n[2]))
    } else {
        return(array(x, dim=n))
    }
}

#' @title Matrix of zeros
#' @description wrapper of `zeros_or_ones()` that replicates the behavior of
#' the `zeros()` function on Matlab
#' @param n1 number of rows
#' @param n2 number of columns
#' @param ... extra dimensions
zeros <- function(n1, n2 = n1, ...) {
    if (length(n1) == 1) {
        n <- c(n1, n2, ...)
    } else {
        n <- n1
    }
    return(zeros_or_ones(n, 0))
}

#' @title Matrix of ones
#' @description wrapper of `zeros_or_ones()` that replicates the behavior of
#' the `ones()` function on Matlab
#' @param n1 number of rows
#' @param n2 number of columns
#' @param ... extra dimensions
ones <- function(n1, n2 = n1, ...) {
    if (length(n1) == 1) {
        n <- c(n1, n2, ...)
    } else {
        n <- n1
    }
    return(zeros_or_ones(n, 1))
}