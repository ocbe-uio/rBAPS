#' @title Convert Matlab function to R
#' @description Performs basic syntax conversion from Matlab to R
#' @param filename name of the file
#' @param saveOutput if `TRUE`, `filename` is overwritten. Defaults to `FALSE`
#' @return text converted to R, printed to screen or replacing input file
#' @author Waldir Leoncio
#' @importFrom utils write.table
#' @export
matlab2r <- function(filename, saveOutput = FALSE) {

	# ======================================================== #
	# Verification                                             #
	# ======================================================== #
	if (!file.exists(filename)) stop("File not found")

	# ======================================================== #
	# Reading file into R                                      #
	# ======================================================== #
	txt <- readLines(filename)

	# ======================================================== #
	# Replacing text                                           #
	# ======================================================== #

	# Function header ---------------------------------------- #
	txt <- gsub(
		pattern     = "function (.+)\\s+=\\s*(.+)\\((.+)\\)",
		replacement = "\\2 <- function(\\3) { return(\\1)",
		x           = txt
	)
	txt <- gsub(
		pattern     = "function (.+)\\((.+)\\)",
		replacement = "\\1 <- function(\\2) {",
		x           = txt
	)
	# txt <- gsub("\\%\\s*(\\w+)", "# \\1", txt)
	txt <- gsub(";", "", txt)
	txt <- gsub("for (.+)=(.+)", "for (\\1 in \\2) {", txt)
	txt <- gsub("end", "}", txt)
	txt <- gsub("(.),(\\S)", "\\1, \\2", txt)
	# TODO: replace forms like (:,:) with [, ] if they come before <-
	# TODO: add argument to skip some of these rules
	txt <- gsub("if (.+)", "if (\\1) {", txt) # FIXME: paste comments after {
	txt <- gsub("else$", "} else {", txt)
	txt <- gsub("elseif", "} else if", txt)
	txt <- gsub("\\(~", "(!", txt)
	txt <- gsub("while (.+)", "while \\1 {", txt)
	## Math operators
	txt <- gsub("(\\S)\\+(\\S)", "\\1 + \\2", txt)
	txt <- gsub("(\\S)\\-(\\S)", "\\1 - \\2", txt)
	txt <- gsub("(\\S)\\*(\\S)", "\\1 * \\2", txt)
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
