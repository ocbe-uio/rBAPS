#' @title Prompt for multiple-choice
#' @param quest Question
#' @param dlgtitle Title of question
#' @param btn Vector of alternatives
#' @param defbtn Scalar with the name of the default option
#' @param accepted_ans Vector containing accepted answers
#' @description This function aims to loosely mimic the behavior of the
#' questdlg function on Matlab
#' @export
questdlg <- function(quest,
                     dlgtitle = "",
                     btn = c("y", "n"),
                     defbtn = "n",
                     accepted_ans = c("y", "yes", "n", "no")) {
  message(dlgtitle)
  # ==========================================================================
  # Replacing the default option with a capitalized version on btn
  # ==========================================================================
  btn[match(tolower(defbtn), tolower(btn))] <- toupper(defbtn)
  # ==========================================================================
  # Creating prompt
  # ==========================================================================
  option_char <- paste0(" [", paste(btn, collapse = ", "), "]")
  answer <- readline(paste0(quest, option_char, ": "))
  # ==========================================================================
  # Processing answer
  # ==========================================================================
  answer <- tolower(answer)
  if (!(answer %in% tolower(c(btn, accepted_ans)))) {
    if (answer != "") {
      warning(
        "'", answer, "' is not a valid alternative. Defaulting to ",
        defbtn
      )
    }
    answer <- defbtn
  }
  return(answer)
}
