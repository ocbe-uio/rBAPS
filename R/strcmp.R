#' @title Compare two character elements
#' @description Logical test if two character elements are identical
#' @param s1 first character element (string, vector or matrix)
#' @param s2 second character element (string, vector or matrix)
#' @return a logical element of the same type as the input
#' @export
strcmp <- function(s1, s2) {
  if (length(s1) == 1 & length(s2) == 1) {
    # Both are scalars, comparison is straightforward
    return(identical(s1, s2))
  } else if (length(s1) == 1 & length(s2) > 1) {
    # s1 is a scalar and s2 is a vector or a matrix
    checks <- sapply(s2, function(s) s1 %in% s)
    if (is(s2, "matrix")) checks <- matrix(checks, nrow(s2))
  } else if (length(s1) > 1 & length(s2) == 1) {
    # s1 is a vector/matrix, s2 is a scalar
    checks <- sapply(s1, function(s) s2 %in% s)
    if (is(s1, "matrix")) checks <- matrix(checks, nrow(s1))
  } else {
    # s1 and s2 are vectors/matrices
    if (identical(dim(s1), dim(s2))) {
      checks <- as.matrix(s1 == s2)
    } else {
      stop("Inputs must be the same size or either one can be a scalar.")
    }
  }
  return(checks)
}
