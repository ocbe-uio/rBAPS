#' @title Cell array
#' @description Creates an array of zeros
#' @param n a the first dimension (or both, if sz is not passed)
#' @param sz the second dimension (or 1st and 2nd, if not passed)
#' @return An array of zeroes with the dimensions passed on call
cell <- function(n, sz = c(n, n), ...) {
	if (length(sz) == 1 & missing(...)) {
		return(array(dim = c(n, sz)))
	} else if (length(sz) == 2) {
		return(array(dim = sz))
	} else {
		return(array(dim = c(n, sz, ...)))
	}
}