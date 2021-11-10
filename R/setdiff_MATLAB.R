#' @title Set differences of two arrays
#' @description Loosely replicates the behavior of the homonym Matlab function
#' @param A first array
#' @param B second array
#' @param legacy if `TRUE`, preserves the behavior of the setdiff function from MATLAB R2012b and prior releases. (currently not supported)
#' @author Waldir Leoncio
setdiff_MATLAB <- function(A, B, legacy = FALSE) {
  if (legacy) message("legacy=TRUE not supported. Ignoring.")
  if (is(A, "numeric") & is(B, "numeric")) {
    values <- sort(unique(A[is.na(match(A, B))]))
  } else if (is(A, "data.frame") & is(B, "data.frame")) {
    C <- A
    exclude_rows <- vector()
    for (r1 in seq_len(nrow(A))) {
      for (r2 in seq_len(nrow(B))) {
        if (all(A[r1, ] == B[r2, ])) {
          exclude_rows <- append(exclude_rows, r1)
        }
      }
    }
    values <- C[-exclude_rows, ]
  }
  # TODO: add support for indices (if necessary)
  return(values)
}
