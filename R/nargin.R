#' @title Number of function input arguments
#' @description Returns the number of arguments passed to the parent function
#' @return An integer
#' @author Waldir Leoncio
#' @note This function only makes sense inside another function
nargin <- function() {
	# FIXME: returning 0 because it is using its own envir instead of parent's
	print(parent.env(environment()))
	length(as.list(match.call(envir = parent.env(environment())))) - 1
	# length(ls(envir=parent.env(environment()))) - 1
}