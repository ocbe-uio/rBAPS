#' @title Simulate All Frequencies
#' @description Lisää jokaista alleelia joka populaation joka lokukseen j1/noalle(j) verran. Näin saatuja counts:eja vastaavista Dirichlet-jakaumista simuloidaan arvot populaatioiden alleelifrekvensseille.
#' Add each allele to each locus in each population by j 1 / noalle(j). The Dirichlet distributions corresponding to the counts thus obtained simulate values for the allele frequencies of the populations.
#' @param noalle noalle
simulateAllFreqs <- function(noalle) {
  if (isGlobalEmpty(globals$COUNTS)) {
    max_noalle <- 0
    nloci <- 0
    npops <- 1
  } else {
    max_noalle <- size(globals$COUNTS, 1)
    nloci <- size(globals$COUNTS, 2)
    npops <- size(globals$COUNTS, 3)
  }

  prioriAlleelit <- zeros(max_noalle, nloci)
  if (nloci > 0) {
    for (j in 1:nloci) {
      prioriAlleelit[1:noalle[j], j] <- 1 / noalle[j]
    }
  }
  prioriAlleelit <- repmat(prioriAlleelit, matrix(c(1, 1, npops), 1))
  counts <- ifelse(
    test = isGlobalEmpty(globals$COUNTS),
    yes  = prioriAlleelit,
    no   = globals$COUNTS + prioriAlleelit
  )
  allfreqs <- zeros(size(counts))

  for (i in 1:npops) {
    if (nloci > 0) {
      for (j in 1:nloci) {
        simuloidut <- randdir(counts[1:noalle[j], j, i], noalle[j])
        allfreqs[1:noalle[j], j, i] <- simuloidut
      }
    }
  }
  return(allfreqs)
}
