mysubset <- function(small, large) {
  #  MYSUBSET Is the small set of + ve integers a subset of the large set?
  #  p <- mysubset(small, large)

  #  Surprisingly, this is not built - in.

  if (is.null(small)) {
    p <- 1# is.null(large)
  } else {
    p <- length(myintersect(small, large)) == length(small)
  }
  return(p)
}
