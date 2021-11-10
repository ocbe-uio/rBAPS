#' @title Blanks
#' @description Create character vector of blanks
#' @details This function emulates the behavior of a homonimous function from Matlab
#' @param n length of vector
#' @return Vector of n blanks
#' @author Waldir Leoncio
#' @export
blanks <- function(n) {
  if (n < 0) {
    warning("Negative n passed. Treating as n = 0")
    n <- 0
  }
  paste(rep(" ", n), collapse = "")
}
