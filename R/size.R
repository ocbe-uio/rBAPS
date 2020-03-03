#' @title Size of an object
#' @description This functions tries to replicate the behavior of the base function "size" in Matlab
#' @param x object to be evaluated
#' @param d dimension of object to be evaluated
#' @note On MATLAB, size(1, 100) returns 1. As a matter of fact, if the user
#' calls for a dimension which x doesn't have `size()` always returns 1. R's
#' default behavior is more reasonable in those cases (i.e., returning NA),
#' but since the point of this function is to replicate MATLAB behaviors
#' (bugs and questionable behaviors included), this function also does this.
#' @export
size <- function(x, d) {
    # Determining the number of dimensions
    if (all(is.na(x))) {
        if (missing(d)) {
            return(c(0, 0))
        } else {
            return(ifelse(d <= 2, 0, 1))
        }
    }
    if (length(x) == 1) {
        # x is surely a scalar
        return(1)
    } else {
        # x is a vector, a matrix or an array
        n_dim <- ifelse(is.null(dim(x)), 1, length(dim(x)))
        if (missing(d)) {
            if (n_dim == 1) {
                out <- c(1, length(x))
            } else {
                out <- dim(x)
            }
        } else {
            out <- ifelse(n_dim == 1, c(1, length(x))[d], dim(x)[d])
            if (is.na(out)) out <- 1  # for MATLAB compatibility
        }
        return(out)
    }
}