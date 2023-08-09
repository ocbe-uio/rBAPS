indMixWrapper <- function(c, npops, fixedK = FALSE) {
  if (fixedK) {
    stop("indMix_fixK() not yet implemented.") # TODO: translate indMix_fixK.m
  } else {
    indMix(c, npops, TRUE)
  }
}
