#' @title Vector creation
#' @description Simulates the function `colon()` and its equivalent `:` operator from Matlab, which have a similar but not quite equivalent behavior when compared to `seq()` and `:` in R.
#' @param a initial number
#' @param b final number
#' @export
colon <- function(a, b) {
  if (a <= b) {
    return(a:b)
  } else {
    return(vector(mode = "numeric"))
  }
}
