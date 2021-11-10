#' @title Find indices and values of nonzero elements
#' @description Emulates behavior of `find`
#' @param x object or logic operation on an object
#' @param sort sort output?
find <- function(x, sort = TRUE) {
  if (is.logical(x)) {
    out <- which(x)
  } else {
    out <- which(x > 0)
  }
  if (sort) {
    out <- sort(out)
  }
  return(out)
}
