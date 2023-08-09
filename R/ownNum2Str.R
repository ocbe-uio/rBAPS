#' @title Own number to string
#' @description Converts numbers to strings
#' @param number number
#' @note On Matlab, if number is NaN the output is 'NaN'. Here, the output will be an error. Also, the function belo expects "number" to have length one, whereas Matlab accepts vectors.
ownNum2Str <- function(number) {
  absolute <- abs(number)
  if (absolute < 1000) {
    str <- as.character(number)
  } else if (absolute < 10000000) {
    first_three <- number %% 1000
    next_four <- (number - first_three) / 1000
    first_three <- abs(first_three)
    if (first_three < 10) {
      first_three <- paste0("00", as.character(first_three))
    } else if (first_three < 100) {
      first_three <- paste0("0", as.character(first_three))
    } else {
      first_three <- as.character(first_three)
    }
    str <- paste0(as.character(next_four), first_three)
  } else if (absolute < 100000000) {
    first_four <- number %% 10000
    next_four <- (number - first_four) / 10000
    first_four <- abs(first_four)
    if (first_four < 10) {
      first_four <- paste0("000", as.character(first_four))
    } else if (first_four < 100) {
      first_four <- paste0("00", as.character(first_four))
    } else if (first_four < 1000) {
      first_four <- paste0("0", as.character(first_four))
    } else {
      first_four <- as.character(first_four)
    }
    str <- paste0(as.character(next_four), first_four)
  } else {
    str <- as.character(number)
  }
  return(str)
}
