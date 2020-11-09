#' @title Set differences of two arrays
#' @description Loosely replicates the behavior of the homonym Matlab function
#' @param A first array
#' @param B second array
#' @param legacy if `TRUE`, preserves the behavior of the setdiff function from MATLAB R2012b and prior releases. (currently not supported)
#' @author Waldir Leoncio
setdiff_MATLAB <- function(A, B, legacy = FALSE) {
	if (legacy) message("legacy=TRUE not supported. Ignoring.")
	if (is(A, "numeric") & is(B, "numeric")) {
		values <- sort(unique(A[is.na(match(A, B))]))
	} else if (is(A, "data.frame") & is(B, "data.frame")) {
		stop("Not implemented for data frames")
	}
	# TODO: add support for indices (if necessary)
	return(values)
}