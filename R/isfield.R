#' @title Checks if a list contains a field
#' @description This function tries to replicate the behavior of the `isfield`
#' function in Matlab
#' @param x list
#' @param field name of field
#' @references https://se.mathworks.com/help/matlab/ref/isfield.html
isfield <- function(x, field) {
    field %in% names(x)
}