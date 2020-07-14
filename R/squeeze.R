#' @title Squeeze
#' @description Remove dimensions of length 1
#' @details This function implements the behavior of the homonimous function on
#' Matlab. `B = squeeze(A)` returns an array with the same elements as the
#' input array A, but with dimensions of length 1 removed. For example, if A is
#' a 3-by-1-by-1-by-2 array, then squeeze(A) returns a 3-by-2 matrix. If A is a
#' row vector, column vector, scalar, or an array with no dimensions of length
#' 1, then squeeze returns the input A.
#' @note This is basically a wrapper of drop() with a minor adjustment to adapt
#' the output to what happens on Matlab
#' @param A input or array matrix
#' @return An array with the same elements as the input array, but with
#' dimensions of length 1 removed.
#' @author Waldir Leoncio
squeeze <- function(A) as.matrix(drop(A))