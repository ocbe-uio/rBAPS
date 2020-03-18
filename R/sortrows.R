#' @title Sort rows of matrix or table
#' @description Emulates the behavior of the `sortrows` function on Matlab
#' @param A matrix
#' @param column ordering column
sortrows <- function(A, column = 1) {
    if (length(column) == 1) {
        new_row_order <- order(A[, column])
    } else if (length(column) == 2) {
        new_row_order <- order(A[, column[1]], A[, column[2]])
    } else {
        stop("Not yet implemented for 2+ tie-breakers")
    }
    A_reordered <- A[new_row_order, ]
    return(A_reordered)
}