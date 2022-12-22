myintersect <- function(A, B) {
  #  MYINTERSECT Intersection of two sets of positive integers (much faster than built - in intersect)
  #  C <- myintersect(A, B)

  A <- t(A)
  B <- t(B)

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

  if (ma == 0 | mb == 0) {
    C <- vector()
  } else {
    # bits <- sparse(1, max(ma, mb))
    bits <- zeros(1, max(ma, mb))
    bits[A] <- 1
    C <- B[as.logical(bits[B])]
  }
  return(C)
}
