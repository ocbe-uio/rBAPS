#' @title Squeeze
#' @description Remove dimensions of length 1
#' @details This function implements the behavior of the homonimous function on
#' Matlab. `B = squeeze(A)` returns an array with the same elements as the
#' input array A, but with dimensions of length 1 removed. For example, if A is
#' a 3-by-1-by-1-by-2 array, then squeeze(A) returns a 3-by-2 matrix. If A is a
#' row vector, column vector, scalar, or an array with no dimensions of length
#' 1, then squeeze returns the input A.
#' @param A input or array matrix
#' @return An array with the same elements as the input array, but with
#' dimensions of length 1 removed.
#' @author Waldir Leoncio
squeeze <- function(A) {
	A <- as.array(A)
	dim_1 <- which(dim(A) == 1)
	B <- array(A, dim = dim(A)[-dim_1])

	# Workaround to match Matlab behavior
	if (length(dim(B)) == 1) B <- as.matrix(B)

	return(B)
}