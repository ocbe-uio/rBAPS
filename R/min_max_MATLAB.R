#' @title Minimum (MATLAB version)
#' @description Finds the minimum value for each column of a matrix, potentially returning the indices instead
#' @param X matrix
#' @param indices return indices?
#' @return Either a list or a vector
#' @author Waldir Leoncio
min_MATLAB <- function(X, indices = TRUE) {
  if (!is(X, "matrix")) X <- as.matrix(X)
  mins <- apply(X, 2, min)
  idx <- sapply(seq_len(ncol(X)), function(x) match(mins[x], X[, x]))
  if (indices) {
    return(list(mins = mins, idx = idx))
  } else {
    return(mins)
  }
}

#' @title Maximum (MATLAB version)
#' @description Finds the minimum value for each column of a matrix, potentially returning the indices instead
#' @param X matrix
#' @param indices return indices?
#' @return Either a list or a vector
#' @author Waldir Leoncio
max_MATLAB <- function(X, indices = TRUE) {
  if (!is(X, "matrix")) X <- as.matrix(X)
  maxs <- apply(X, 2, max)
  idx <- sapply(seq_len(ncol(X)), function(x) match(maxs[x], X[, x]))
  if (indices) {
    return(list(maxs = maxs, idx = idx))
  } else {
    return(maxs)
  }
}
