#' @title Read line from file, removing newline characters
#' @description Equivalent function to its homonymous Matlab equivalent.
#' @param file character vector to be read, usually an output of `fopen()`
#' @return If the file is nonempty, then fgetl returns tline as a character vector. If the file is empty and contains only the end-of-file marker, then fgetl returns tline as a numeric value -1.
#' @author Waldir Leoncio
#' @seealso fopen
#' @export
fgetl <- function(file) {
	# ==========================================================================
	# Validation
	# ==========================================================================
	if (length(file) <= 1) return(-1)
	# ==========================================================================
	# Returning file minus the first line
	# ==========================================================================
	out <- file[-1]
	return(out)
}

#' @title Open file
#' @description Open a text file
#' @param filename Path and name of file to be open
#' @return The same as `readLines(filename)`
#' @author Waldir Leoncio
#' @seealso fgetl
#' @export
fopen <- function(filename) readLines(filename)