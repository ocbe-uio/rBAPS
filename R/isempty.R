#' @title Is Array Empty?
#' @description Determine whether array is empty. An empty array, table, or timetable has at least one dimension with length 0, such as 0-by-0 or 0-by-5.
#' @details Emulates the behavior of the `isempty` function on Matlab
#' @param x array
#'
isempty <- function(x) {
  if (class(x)[1] %in% c("array", "matrix")) {
    dim_mat_x <- dim(x)
  } else {
    dim_mat_x <- dim(matrix(x))
  }
  return(any(dim_mat_x == 0) | is.null(dim_mat_x))
}
