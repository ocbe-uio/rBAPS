#' @title Convert proportion to string
#' @param prob belongs to [0.00, 0.01, ... ,1]
#' @return a 4-mark presentation of proportion
#' @note The `round` function in R, being ISO-compliant, rounds 8.5 to 8. The
#' Matlab equivalent rounds it to 9.
#' @export
proportion2str <- function(prob) {
  if (abs(prob) < 1e-3) {
    str <- "0.00"
  } else if (abs(prob - 1) < 1e-3) {
    str <- "1.00"
  } else {
    prob <- round(100 * prob)
    if (prob < 10) {
      str <- paste0("0.0", as.character(prob))
    } else {
      str <- paste0("0.", as.character(prob))
    }
  }
  return(str)
}
