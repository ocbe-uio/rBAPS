#' @title Element-wise matrix multiplication
#' @description Emulates the `times()` and `.*` operators from Matlab.
#' @details This function basically handles elements of different length better than the `*` operator in R, at least as far as behavior from a Matlab user is expecting.
#' @param a first factor of the multiplication
#' @param b second factor of the multiplication
#' @export
#' @returns matrix with dimensions equal to the larger of the two factors
times <- function(a, b) {
    # Converting everything to matrix because Matlab looooooves the matrix
    a <- as.matrix(a)
    b <- as.matrix(b)

    dominant_mx <- NULL
    if (!all(dim(a) == dim(b))) {
        if (all(dim(a) >= dim(b))) {
            dominant_mx <- a
            dominated_mx <- b
        }  else if (all(dim(b) >= dim(a))) {
            dominant_mx <- b
            dominated_mx <- a
        } else {
            dominant_mx <- "neither"
            dominant_dim <- c(max(nrow(b), nrow(a)), max(ncol(b), ncol(a)))
        }
    }

    if (is.null(dominant_mx)) {
        out <- a * b
    } else if (dominant_mx[1] == "neither") {
        a <- repmat(
            mx = a,
            n = c(dominant_dim[1] - nrow(a) + 1, dominant_dim[2] - ncol(a) + 1)
        )
        b <- repmat(
            mx = b,
            n = c(dominant_dim[1] - nrow(b) + 1, dominant_dim[2] - ncol(b) + 1)
        )
        out <- a * b
    } else {
        # Expanding dominated matrix
        dominated_mx <- repmat(
            mx = dominated_mx,
            n = c(
                nrow(dominant_mx) - nrow(dominated_mx) + 1,
                ncol(dominant_mx) - ncol(dominated_mx) + 1
            )
        )
        out <- dominant_mx * dominated_mx
    }
    return(out)
}