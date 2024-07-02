returnInOrder <- function(inds, pop, globalRows, data, adjprior, priorTerm) {
  # % Palauttaa yksil�t j�rjestyksess� siten, ett� ensimm�isen� on
  # % se, jonka poistaminen populaatiosta pop nostaisi logml:n
  # % arvoa eniten.

  ninds <- length(inds)
  apuTaulu <- cbind(inds, zeros(ninds, 1))

  for (i in 1:ninds) {
    ind <- inds[i]
    rows <- globalRows[i, 1]:globalRows[i, 2]
    diffInCounts <- computeDiffInCounts(
      rows, size(globals$COUNTS, 1), size(globals$COUNTS, 2), data
    )
    diffInSumCounts <- colSums(diffInCounts)

    globals$COUNTS[, , pop] <- globals$COUNTS[, , pop] - diffInCounts
    globals$SUMCOUNTS[pop, ] <- globals$SUMCOUNTS[pop, ] - diffInSumCounts
    apuTaulu[i, 2] <- computePopulationLogml(pop, adjprior, priorTerm)
    globals$COUNTS[, , pop] <- globals$COUNTS[, , pop] + diffInCounts
    globals$SUMCOUNTS[pop, ] <- globals$SUMCOUNTS[pop, ] + diffInSumCounts
  }
  apuTaulu <- sortrows(apuTaulu, 2)
  inds <- apuTaulu[ninds:1, 1]
  return(inds)
}
