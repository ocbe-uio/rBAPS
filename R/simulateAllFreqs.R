#' @title Simulate All Frequencies
#' @description Lisää jokaista alleelia joka populaation joka lokukseen j1/noalle(j) verran. Näin saatuja counts:eja vastaavista Dirichlet-jakaumista simuloidaan arvot populaatioiden alleelifrekvensseille.
#' Add each allele to each locus in each population by j 1 / noalle(j). The Dirichlet distributions corresponding to the counts thus obtained simulate values for the allele frequencies of the populations.
#' @param noalle noalle
#' @param COUNTS COUNTS
#' @export

simulateAllFreqs <- function(noalle, COUNTS = matrix(NA, 0, 0)) {
    max_noalle <- size(COUNTS, 1)
    nloci <- size(COUNTS, 2)
    npops <- size(COUNTS, 3)

    prioriAlleelit <- zeros(max_noalle, nloci)
    if (nloci > 0) {
        for (j in 1:nloci) {
            prioriAlleelit[1:noalle[j], j] <- 1 / noalle[j]
        }
    }
    prioriAlleelit <- repmat(prioriAlleelit, matrix(c(1, 1, npops), 1))
    counts <- COUNTS + prioriAlleelit
    allfreqs <- zeros(size(counts))

    for (i in 1:npops) {
        if (nloci > 0) {
            for (j in 1:nloci) {
                simuloidut <- randdir(counts[1:noalle[j], j, i] , noalle[j])
                allfreqs[1:noalle[j], j, i] <- simuloidut
            }
        }
    }
    return(allfreqs)
}