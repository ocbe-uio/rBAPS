#' @title Prompt for multiple-choice
#' @param quest Question
#' @param dlgtitle Title of question
#' @param btn Vector of alternatives
#' @param defbtn Scalar with the name of the default option
#' @description This function aims to loosely mimic the behavior of the
#' questdlg function on Matlab
#' @export
questdlg <- function(quest, dlgtitle, btn = c('y', 'n'), defbtn = 'n') {
	message(dlgtitle)
	# ==========================================================================
	# Replacing the default option with a capitalized version on btn
	# ==========================================================================
	btn[match(tolower(defbtn), tolower(btn))] <- toupper(defbtn)
	# ==========================================================================
	# Creating prompt
	# ==========================================================================
	option_char <- paste0(' [', paste(btn, collapse = ', '), ']')
	answer <- readline(paste0(quest, option_char, ": "))
	# ==========================================================================
	# Processing answer
	# ==========================================================================
	answer <- tolower(answer)
	if (!(answer %in% tolower(c(btn)))) {
		if (answer != "") {
			warning(
				"'", answer, "' is not a valid altenative. Defaulting to ",
				defbtn
			)
		}
		answer <- defbtn
	}
	return(answer)
}