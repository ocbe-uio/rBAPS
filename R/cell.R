#' @title Cell array
#' @description Creates an array of zeros
#' @param n a the first dimension (or both, if sz is not passed)
#' @param sz the second dimension (or 1st and 2nd, if not passed)
#' @param expandable if TRUE, output is a list (so it can take different
#' lengths)
#' @param ... Other dimensions
#' @return An array of zeroes with the dimensions passed on call
cell <- function(n, sz = c(n, n), expandable=FALSE, ...) {

	# Uglyly figuring out if the third arg is an extra dim --- #

	sz3 <- vector()
	if (!is.logical(expandable)) {
		sz3 <- expandable
		expandable <- FALSE
	}
	args <- c(as.list(environment()), list(...))
	exp <- args$expandable
	extra_dims <- c(sz3, args[names(args) == ""])

	# Creating output vector --------------------------------- #

	if (exp) {
		return(vector("list", length = n))
	}
	if (length(sz) == 1 & length(extra_dims) == 0) {
		return(array(0, dim = c(n, sz)))
	} else if (length(extra_dims) > 0) {
		return(array(0, dim = c(n, sz, extra_dims)))
	} else if (length(sz) == 2) {
		return(array(0, dim = sz))
	} else {
		return(array(0, dim = c(n, sz, ...)))
	}
}