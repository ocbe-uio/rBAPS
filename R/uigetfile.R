#' @title Select a file for loading
#' @description Loosely mimics the functionality of the `uigetfile` function on
#' Matlab.
#' @references https://se.mathworks.com/help/matlab/ref/uigetfile.html
#' @param title Pre-prompt message
#' @export
uigetfile <- function(title =  "") {
	# ==========================================================================
	# Pre-prompt message
	# ==========================================================================
	cat(title)
	# ==========================================================================
	# Reading file path and name
	# ==========================================================================
	filepath <- readline(
		paste0("Enter file path (leave empty for ", getwd(), "): ")
	)
	if (filepath == "") filepath <- getwd()
	filename <- file.choose()
	# ==========================================================================
	# Organizing output
	# ==========================================================================
	out <- list(name = filename, path = filepath)
	return(out)
}