#' @title Save file
#' @param filter accepted file extension
#' @param title Title
#' @description This function intends to loosely mimic the behaviour of the
#' homonymous Matlab function.
#' @export
uiputfile <- function(filter = ".rda", title = "Save file") {
	# ==========================================================================
	# Processing input
	# ==========================================================================
	message(title)
	filename <- readline(paste0('File name (end with ', filter, '): '))
	filepath <- readline(paste0('File path (leave empty for ', getwd(), '): '))
	if (filename == "") filename <- 0
	if (filepath == "") filepath <- getwd()
	# ==========================================================================
	# Processing output
	# ==========================================================================
	out <- list(name = filename, path = filepath)
	return(out)
}