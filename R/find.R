#' @title Find indices and values of nonzero elements
#' @description Emulates behavior of `find`
find <- function(x) {
    if (is.logical(x)) {
        return(which(x))
    } else {
        return(which(x > 0))        
    }
}