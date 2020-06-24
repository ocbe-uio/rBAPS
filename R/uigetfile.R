#' @title Select a file for loading
#' @description Loosely mimics the functionality of the `uigetfile` function on
#' Matlab.
#' @references https://se.mathworks.com/help/matlab/ref/uigetfile.html
#' @param filter Filter listed files
#' @param title Pre-prompt message
#' @export
uigetfile <- function(filter = "", title =  "") {
	# ==========================================================================
	# Pre-prompt message
	# ==========================================================================
	message(title)
	# ==========================================================================
	# Reading file path and name
	# ==========================================================================
	filepath <- readline(
		paste0("Enter file path (leave empty for ", getwd(), "): ")
	)
	if (filepath == "") filepath <- getwd()
	# ==========================================================================
	# Presenting possible files
	# ==========================================================================
	message("Files present on that directory:")
	print(list.files(path = filepath, pattern = filter, ignore.case = TRUE))
	filename <- file.choose()
	# ==========================================================================
	# Organizing output
	# ==========================================================================
	out <- list(name = filename, path = filepath)
	return(out)
}