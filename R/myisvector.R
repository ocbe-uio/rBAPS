myisvector <- function(V) {
  # Kuten isvector(V)

  A <- size(V)
  r <- (length(A) == 2) & (min(A) == 1)
  return(r)
}
