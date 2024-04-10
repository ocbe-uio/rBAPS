#' @title Compute all freqs - version 2
#' @description Lisää a priori jokaista alleelia joka populaation joka lokukseen
#' j 1/noalle(j) verran.
#' @param noalle noalle
computeAllFreqs2 <- function(noalle) {
  globals$COUNTS <- ifelse(isGlobalEmpty(globals$COUNTS), vector(), globals$COUNTS)
  globals$COUNTS <- ifelse(isGlobalEmpty(globals$COUNTS), vector(), globals$COUNTS)
  max_noalle <- size(globals$COUNTS, 1)
  nloci <- size(globals$COUNTS, 2)
  npops <- size(globals$COUNTS, 3)
  sumCounts <- globals$COUNTS + ones(size(globals$COUNTS))
  sumCounts <- reshape(t(sumCounts), c(1, nloci, npops))
  sumCounts <- repmat(sumCounts, c(max_noalle, 1, 1))

  prioriAlleelit <- zeros(max_noalle, nloci)
  if (nloci > 0) {
    for (j in 1:nloci) {
      prioriAlleelit[1:noalle[j], j] <- 1 / noalle[j]
    }
  }
  prioriAlleelit <- repmat(prioriAlleelit, c(1, 1, npops))
  counts <- ifelse(
    test = isGlobalEmpty(globals$COUNTS),
    yes  = prioriAlleelit,
    no   = globals$COUNTS + prioriAlleelit
  )
  allFreqs <- counts / drop(sumCounts)
  return(allFreqs)
}
