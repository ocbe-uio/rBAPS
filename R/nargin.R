#' @title Number of function input arguments
#' @description Returns the number of arguments passed to the parent function
#' @return An integer
#' @author Waldir Leoncio
#' @note This function only makes sense inside another function
#' @references https://stackoverflow.com/q/64422780/1169233
nargin <- function() {
  if(sys.nframe() < 2) stop("must be called from inside a function")
  length(as.list(sys.call(-1))) - 1
}
