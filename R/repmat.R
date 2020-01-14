#' @title Repeat matrix
#' @description Repeats a matrix over n columns and rows
#' @details This function was created to replicate the behavior of a homonymous
#' function on Matlab
#' @param mx matrix
#' @param n either a scalar with the number of replications in both rows and columns or a 2-length vector with individual repetitions.
#' @return matrix replicated over `ncol(mx) * n` columns and `nrow(mx) * n` rows
#' @note The Matlab implementation of this function accepts `n` with length > 2.
#' @export
repmat <- function (mx, n) {
    if (length(n) > 2) warning("Extra dimensions of n ignored")
    if (length(n) == 1) n <- rep(n, 2)
    out <- mx_cols <- rep(mx, n[1])
    if (n[2] > 1) {
        for (i in seq(n[2] - 1)) out <- rbind(out, mx_cols)
    }
    return(unname(as.matrix(out)))
}