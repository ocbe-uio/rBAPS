#' @title Read line from file, removing newline characters
#' @description Equivalent function to its homonymous Matlab equivalent.
#' @param file file to be read
#' @return If the file is nonempty, then fgetl returns tline as a character vector. If the file is empty and contains only the end-of-file marker, then fgetl returns tline as a numeric value -1.
#' @author Waldir Leoncio
#' @export
fgetl <- function(file) {
	# ==========================================================================
	# Validation
	# ==========================================================================
	if (file == "") return(-1)
	# ==========================================================================
	# Determine next line to be read
	# ==========================================================================
	if (is.null(attr(file, "last_read_line"))) {
		attr(file, "last_read_line") <- 1
	} else {
		attr(file, "last_read_line") <- attr(file, "last_read_line") + 1
	}
	# ==========================================================================
	# Returning next line
	# ==========================================================================
	out <- file[attr(file, "last_read_line")]
	return(out)
}

#' @title Open file
#' @description Open a text file
#' @param filename Path and name of file to be open
#' @return The same as `readLines(filename)`
#' @author Waldir Leoncio
#' @export
fopen <- function(filename) readLines(filename)