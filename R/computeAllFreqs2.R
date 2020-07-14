#' @title Compute all freqs - version 2
#' @description Lisää a priori jokaista alleelia joka populaation joka lokukseen
#' j 1/noalle(j) verran.
#' @param noalle noalle
#' @param COUNTS counts
#' @param SUMCOUNTS sumcounts
#' @export
computeAllFreqs2 <- function (noalle, COUNTS = matrix(NA, 0, 0),
    SUMCOUNTS = NA) {

    max_noalle <- size(COUNTS, 1)
    nloci <- size(COUNTS,2)
    npops <- size(COUNTS,3)

    sumCounts <- SUMCOUNTS + ones(size(SUMCOUNTS))
    sumCounts <- reshape(t(sumCounts), c(1, nloci, npops))
    sumCounts <- repmat(sumCounts, c(max_noalle, 1, 1))

    prioriAlleelit <- zeros(max_noalle, nloci)
    if (nloci > 0) {
        for (j in 1:nloci) {
            prioriAlleelit[1:noalle[j], j] <- 1 / noalle[j]
        }
    }
    prioriAlleelit <- repmat(prioriAlleelit, c(1, 1, npops))
    counts <- COUNTS + prioriAlleelit
    allFreqs <- counts / drop(sumCounts)
    return(allFreqs)
}