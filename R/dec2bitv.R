dec2bitv <- function(d, n) {
  # DEC2BITV Convert a decimal integer to a bit vector.
  # bits <- dec2bitv(d, n) is just like the built - in dec2bin, except the answer is a vector, not a as.character.
  # n is an optional minimum length on the bit vector.
  # If d is a vector,  each row of the output array will be a bit vector.

  if (!exists("n") || n == 0) {
    n <- 1 # Need at least one digit even for 0.
  }

  f <- e <- NA
  c(f, e) %<-% matlab2r::log2(base::max(d)) # How many digits do we need to represent the numbers?
  bits <- floor(d * 2 ^ (seq(1 - base::max(n, e), 0))) %% 2
  return(bits)
}
