#' @title Set differences of two arrays
#' @description Loosely replicates the behavior of the homonym Matlab function
#' @param A first array
#' @param B second awway
#' @param legacy if `TRUE`, preserves the behavior of
#' @return
#' @author Waldir Leoncio
#' @export
setdiff <- function(A, B, legacy = FALSE) {
	values <- sort(unique(A[is.na(match(A, B))]))
	# browser() # TEMP
	return(values)
}