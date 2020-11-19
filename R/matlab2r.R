#' @title Convert Matlab function to R
#' @description Performs basic syntax conversion from Matlab to R
#' @param filename name of the file
#' @param output can be "asis", "clean" (default) or "save"
#' @param improve_formatting if `TRUE` (default), makes minor changes
#' to conform to best-practice formatting conventions
#' @return text converted to R, printed to screen or replacing input file
#' @author Waldir Leoncio
#' @importFrom utils write.table
#' @export
matlab2r <- function(
	filename, output = "clean", improve_formatting=TRUE
) {
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

	# Uncommenting ------------------------------------------- #
	txt <- gsub("^#\\s?(.+)", "\\1", txt)

	# Function header ---------------------------------------- #
	out <- gsub(
		pattern     = "\\t*function (.+)\\s*=\\s*(.+)\\((.+)\\)",
		replacement = "\treturn(\\1)",
		x           = txt[1]
	)
	txt <- gsub(
		pattern     = "\\t*function (.+)\\s*=\\s*(.+)\\((.+)\\)",
		replacement = "\\2 <- function(\\3) {",
		x           = txt
	)
	txt <- gsub(
		pattern     = "function (.+)\\((.+)\\)",
		replacement = "\\1 <- function(\\2) {",
		x           = txt
	)

	# Function body ------------------------------------------ #
	txt <- gsub("(.+)\\.\\.\\.", "\\1", txt)
	txt <- gsub(";", "", txt)

	# Loops and if-statements
	txt <- gsub("for (.+)=(.+)", "for (\\1 in \\2) {", txt)
	txt <- gsub("end$", "}", txt)
	txt <- gsub("if (.+)", "if (\\1) {", txt) # FIXME: paste comments after {
	txt <- gsub("else$", "} else {", txt)
	txt <- gsub("elseif", "} else if", txt)
	txt <- gsub("while (.+)", "while \\1 {", txt)

	# MATLAB-equivalent functions in R
	txt <- gsub("gamma_ln", "log_gamma", txt)

	# Subsets ------------------------------------------------ #
	txt <- gsub("([^\\(]+)\\((.+)\\)\\s?=(.+)", "\\1[\\2] <- \\3", txt)

	# Formatting --------------------------------------------- #
	if (improve_formatting) {
		txt <- gsub("(.),(\\S)", "\\1, \\2", txt)
		# Math operators
		txt <- gsub("(\\S)\\+(\\S)", "\\1 + \\2", txt)
		txt <- gsub("(\\S)\\-(\\S)", "\\1 - \\2", txt)
		txt <- gsub("(\\S)\\*(\\S)", "\\1 * \\2", txt)
		# Logic operators
		txt <- gsub("\\(~", "(!", txt)
		# Assignment
		txt <- gsub("(.+)\\s?=\\s?(.+)", "\\1 <- \\2", txt)
	}

	# Adding output and end-of-file brace -------------------- #
	txt <- append(txt, paste(out, "\n}"))

	# Returning converted code ------------------------------- #
	if (output == "asis") {
		return(txt)
	} else if (output == "clean") {
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
