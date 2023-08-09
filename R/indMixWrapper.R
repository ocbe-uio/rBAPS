indMixWrapper <- function(fixedK = FALSE) {
  if (fixedK) {
    indMix(c, npops, TRUE)
  } else {
    stop("indMix_fixK() not yet implemented.") # TODO: translate indMix_fixK.m
  }
}
