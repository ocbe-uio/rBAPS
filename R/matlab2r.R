#' @title Convert Matlab function to R
#' @description Performs basic syntax conversion from Matlab to R
#' @param filename name of the file
#' @param output can be "asis", "clean" (default), "save" or "append"
#' @param improve_formatting if `TRUE` (default), makes minor changes
#' to conform to best-practice formatting conventions
#' @param change_assignment if `TRUE` (default), uses `<-` as the assignment operator
#' @param append if `FALSE` (default), overwrites file; otherwise, append
#' output to input
#' @return text converted to R, printed to screen or replacing input file
#' @author Waldir Leoncio
#' @importFrom utils write.table
#' @export
#' @note This function is intended to expedite the process of converting a
#' Matlab function to R by making common replacements. It does not have the
#' immediate goal of outputting a ready-to-use function. In other words,
#' after using this function you should go back to it and make minor changes.
#'
#' It is also advised to do a dry-run with `output = "clean"` and only switching
#' to `output = "save"` when you are confident that no important code will be
#' lost (for shorter functions, a careful visual inspection should suffice).
matlab2r <- function(
	filename, output = "clean", improve_formatting=TRUE, change_assignment=TRUE,
	append=FALSE
) {
	# TODO: this function is too long! Split into subfunctions
	# (say, by rule and/or section)
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
		pattern     = "\\t*function (\\S+)\\s*=\\s*(.+)\\((.+)\\)",
		replacement = "\treturn(\\1)",
		x           = txt[1]
	) # TODO: improve by detecting listed outputs
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
	txt <- gsub("nchoosek", "choose", txt)
	txt <- gsub("isempty", "is.null", txt)
	# txt <- gsub("(.+)\\'", "t(\\1)", txt)

	# Subsets ------------------------------------------------ #
	ass_op <- ifelse(change_assignment, "<-", "=")
	txt <- gsub(
		pattern = "([^\\(]+)\\(([^\\(]+)\\)=(.+)",
		replacement = paste0("\\1[\\2] ", ass_op, "\\3"),
		x = txt
	)
	txt <- gsub("\\(:\\)", "[, ]", txt)
	txt <- gsub("(.+)(\\[|\\():,end(\\]|\\()", "\\1[, ncol()]", txt)

	# Formatting --------------------------------------------- #
	if (improve_formatting) {
		txt <- gsub("(.),(\\S)", "\\1, \\2", txt)
		# Math operators
		txt <- gsub("(\\S)\\+(\\S)", "\\1 + \\2", txt)
		txt <- gsub("(\\S)\\-(\\S)", "\\1 - \\2", txt)
		txt <- gsub("(\\S)\\*(\\S)", "\\1 * \\2", txt)
		# Logic operators
		txt <- gsub("~", "!", txt)
		txt <- gsub("(\\S)>=(\\S)", "\\1 >= \\2", txt)
		txt <- gsub("(\\S)<=(\\S)", "\\1 <= \\2", txt)
		txt <- gsub("(\\S)==(\\S)", "\\1 == \\2", txt)
		# Assignment
		txt <- gsub(
			pattern = "(\\w)(\\s?)=(\\s?)(\\w)",
			replacement = paste0("\\1 ", ass_op, " \\4"),
			x = txt
		)
		# txt <- gsub(
		# 	pattern = "(\\s+(.|\\_|\\[|\\])+)(\\s?)=(\\s?)(.+)",
		# 	replacement = paste0("\\1 ", ass_op, "\\5"),
		# 	x = txt
		# )
	}

	# Adding output and end-of-file brace -------------------- #
	txt <- append(txt, paste(out, "\n}"))

	# Returning converted code ------------------------------- #
	if (output == "asis") {
		return(txt)
	} else if (output == "clean") {
		return(cat(txt, sep="\n"))
	} else if (output == "save") {
		return(
			write.table(
				x         = txt,
				file      = filename,
				quote     = FALSE,
				row.names = FALSE,
				col.names = FALSE,
				append    = append
			)
		)
	} else {
		stop ("Invalid output argument")
	}
}
