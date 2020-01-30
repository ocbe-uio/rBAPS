#' @title Compute rows
#' @description Individuals inds have been given. The function returns a vector,
#' containing the indices of the rows, which contain data from the individuals.
#' @param rowsFromInd rowsFromInd
#' @param inds matrix
#' @param ninds ninds
#' @export
computeRows <- function(rowsFromInd, inds, ninds) {
    if (class(inds) != "matrix") inds <- as.matrix(inds)
    if (identical(dim(inds), c(nrow(inds), 1L))) {
        # Special treatment for vectors because R has col vectors by default,
        # whereas Matlab has row vectors by default.
        inds <- t(inds)
        if (ninds == 0) return(matrix(, 1, 0))
    }
    rows <- inds[, rep(1, rowsFromInd)]
    rows <- rows * rowsFromInd
    miinus <- repmat(t((rowsFromInd - 1):0), c(ninds, 1))
    rows <- rows - miinus
    rows <- matrix(t(rows), c(1, rowsFromInd * ninds))
    return(t(rows))
}

