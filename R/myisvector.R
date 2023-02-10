myisvector <- function(V) {
  # Kuten isvector(V)
  V <- as.matrix(V)
  A <- c(nrow(V), ncol(V))

  r <- (base::max(size(A)) == 2) & (base::min(A) == 1)
  return(r)
}
