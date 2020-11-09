#' @title Convert Matlab function to R
#' @description Performs basic syntax conversion
#' @param filename name of the file
#' @param saveOutput if `TRUE`, `filename` is overwritten. Defaults to `FALSE`
#' @return text converted to R
#' @author Waldir Leoncio
#' @export
matlab2r <- function(filename, saveOutput = FALSE) {
	# Verification
	if (!file.exists(filename)) stop("File not found")
	# Reading file into R
	txt <- readLines(filename)
	# Replacing text
	txt <- gsub(";", "", txt)
	txt <- gsub("for (.+)=(.+)", "for (\\1 in \\2) {", txt)
	txt <- gsub("end", "}", txt)
	# Returning converted code
	if (!saveOutput) {
		return(cat(txt, sep="\n"))
	} else {
		return(
			write.table(
				x         = txt,
				file      = filename,
				quote     = FALSE,
				row.names = FALSE,
				col.names = FALSE
			)
		)
	}
}