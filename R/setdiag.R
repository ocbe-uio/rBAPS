setdiag <- function(M, v) {
  #  SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
  #  M <- set_diag(M, v)

  n <- length(M)
  if (length(v) == 1) {
    v <- repmat(v, c(1, n))
  }

  #  e.g., for 3x3 matrix,  elements are numbered
  #  1 4 7
  #  2 5 8
  #  3 6 9
  #  so diagnoal = [1 5 9]

  J <- seq(1, n ^ 2, n + 1)
  M[J] <- v
  return(M)
}
