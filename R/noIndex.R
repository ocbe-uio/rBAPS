#' @title No index
#' @description Checks that the data contains no index column.
#' @details As input, this function takes two variables from a mixture/admixture
#' result structure.
#' @return puredata: a data contains no index column.
#' @param data data
#' @param noalle noalle
noIndex <- function(data, noalle) {
  limit <- ifelse(is(noalle, "matrix"), ncol(noalle), length(noalle))
  if (size(data, 2) == limit + 1) {
    if (is(data, "matrix")) {
      puredata <- data[, -ncol(data)] # remove the index column
    } else {
      puredata <- data[-length(data)]
    }
  } else {
    puredata <- data
  }
  return(puredata)
}
