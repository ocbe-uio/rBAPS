checkLogml <- function(priorTerm, adjprior, cliques, separators) {
  # tarkistaa logml:n

  # global CLIQCOUNTS
  # global SEPCOUNTS
  # global PARTITION

  npops <- length(unique(PARTITION))
  cliqcounts <- computeCounts(cliques, separators, npops)$cliqcounts
  sepcounts <- computeCounts(cliques, separators, npops)$sepcounts

  CLIQCOUNTS <- cliqcounts
  SEPCOUNTS <- sepcounts

  logml <- computeLogml(adjprior, priorTerm)$logml
  spatialPrior <- computeLogml(adjprior, priorTerm)$spatialPrior

  disp(
    c(
      'logml: ',
      logml2String(logml),
      ', spatial prior: ',
      logml2String(spatialPrior)
    )
  )
  return(logml)
}
