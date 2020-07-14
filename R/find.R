#' @title Find indices and values of nonzero elements
#' @description Emulates behavior of `find`
#' @param x object or logic operation on an object
find <- function(x) {
	if (is.logical(x)) {
		return(which(x))
	} else {
		return(which(x > 0))
	}
}