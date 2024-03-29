myunion <- function(A, B) {
  #  MYUNION Union of two sets of positive integers (much faster than built - in union)
  #  C <- myunion(A, B)

  if (is.null(A)) {
    ma <- 0
  } else {
    ma <- max(A)
  }

  if (is.null(B)) {
    mb <- 0
  } else {
    mb <- max(B)
  }

  if (ma == 0 & mb == 0) {
    C <- vector()
  } else if (ma == 0 & mb > 0) {
    C <- B
  } else if (ma > 0 & mb == 0) {
    C <- A
  } else {
    # bits <- sparse(1, max(ma, mb))
    bits <- zeros(1, max(c(ma, mb)))
    bits[A] <- 1
    bits[B] <- 1
    C <- find(bits)
  }
  return(C)
}
