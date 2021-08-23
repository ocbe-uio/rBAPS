#' @title Check if global variable is empty
#' @description Checks if a global variable has been filled with values other than their initial ones.
#' @details For a list of global variables, check the \code{globals.R} file.
#' @param g the global variable in quesiton.
#' @return TRUE if the variable still contains its original values, FALSE otherwise.
#' @importFrom stats sd
#' @author Waldir Leoncio
isGlobalEmpty <- function(g) {
	return(sum(g) == 0 & sd(g) == 0)
}