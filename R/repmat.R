#' @title Repeat matrix
#' @description Repeats a matrix over n columns and rows
#' @details This function was created to replicate the behavior of a homonymous
#' function on Matlab
#' @param mx matrix
#' @param n either a scalar with the number of replications in both rows and
#' columns or a <= 3-length vector with individual repetitions.
#' @return matrix replicated over `ncol(mx) * n` columns and `nrow(mx) * n` rows
#' @note The Matlab implementation of this function accepts `n` with length > 2.
#'
#' It should also be noted that a concatenated vector in R, e.g. `c(5, 2)`, becomes a column vector when coerced to matrix, even though it may look like a row vector at first glance. This is important to keep in mind when considering the expected output of this function. Vectors in R make sense to be seen as column vectors, given R's Statistics-oriented paradigm where variables are usually disposed as columns in a dataset.
#' @export
repmat <- function(mx, n) {
  # Validation
  if (length(n) > 3) warning("Extra dimensions of n ignored")
  if (!is(mx, "matrix")) mx <- t(as.matrix(mx))
  if (length(n) == 1) n <- rep(n, 2)
  if (any(n == 0)) {
    n_zero <- which(n == 0)
    out_dim <- dim(mx)
    out_dim[n_zero] <- 0
    return(array(dim = out_dim))
  }

  # Replicating cols
  out <- mx_col <- matrix(rep(mx, n[2]), nrow(mx))

  # Replicating rows
  if (n[1] > 1) {
    for (i in seq(n[1] - 1)) out <- rbind(out, mx_col)
  }

  # Replicating 3rd dimension
  if (!is.na(n[3]) & n[3] > 1) out <- array(out, c(dim(out), n[3]))

  # Output
  return(unname(as.array(out)))
}
