indMixWrapper <- function(c, npops, counts, sumcounts, max_iter, fixedK = FALSE, verbose = FALSE) {
  if (fixedK) {
    stop("indMix_fixK() not yet implemented.") # TODO: translate indMix_fixK.m
  } else {
    indMix(c, npops, counts, sumcounts, max_iter, verbose)
  }
}
