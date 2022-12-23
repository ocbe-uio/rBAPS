ind2subv <- function(siz, ndx) {
  #  IND2SUBV Like the built - in ind2sub, but returns the answer as a row vector.
  #  sub <- ind2subv(siz, ndx)
  #  siz and ndx can be row or column vectors.
  #  sub will be of size length(ndx) * length(siz).
  #  Example
  #  ind2subv([2 2 2], 1:8) returns
  #   [1 1 1
  #    2 1 1
  #     ...
  #    2 2 2]
  #  That is, the leftmost digit toggle fastest.
  #
  #  See also SUBV2IND

  n <- length(siz)

  if (n == 0) {
    sub <- ndx
    return(sub)
  }

  if (all(siz == 2)) {
    sub <- dec2bitv(ndx - 1, n)
    sub <- sub[, seq(n, 1, - 1)] + 1
    return(sub)
  }

  cp <- c(1, cumprod(t(siz[])))
  ndx <- ndx[] - 1
  sub <- zeros(length(ndx), n)
  for (i in seq(n, 1, -1)) {# i'th digit
    sub[, i] <- floor(ndx / cp[i]) + 1
    ndx <- ndx %% cp(i)
  }
  return(sub)
}
