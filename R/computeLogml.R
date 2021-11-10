computeLogml <- function(counts, sumcounts, noalle, data, rowsFromInd) {
  nloci <- size(counts, 2)
  npops <- size(counts, 3)
  adjnoalle <- zeros(max(noalle), nloci)
  for (j in 1:nloci) {
    adjnoalle[1:noalle[j], j] <- noalle(j)
    if ((noalle(j) < max(noalle))) {
      adjnoalle[noalle[j] + 1:ncol(adjnoalle), j] <- 1
    }
  }

  rowsInG <- size(data, 1) + rowsFromInd

  logml <- sum(
    sum(
      sum(
        GAMMA_LN[
          counts + 1 +
            repmat(rowsInG * (adjnoalle - 1), c(1, 1, npops))
        ]
      )
    )
  ) -
    npops * sum(sum(GAMMA_LN[1, adjnoalle])) -
    sum(sum(GAMMA_LN[sumcounts + 1, 1]))
  return(logml)
}
