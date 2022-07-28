computeCounts <- function(cliques, separators, npops, PARTITION) {
  ncliq <- size(cliques, 1)
  nsep <- size(separators, 1)

  cliqPartition <- zeros(ncliq, size(cliques, 2))
  sepPartition <- zeros(nsep, size(separators, 2))

  apuCliq <- find(cliques > 0)
  apuSep <- find(separators > 0)

  cliqPartition[apuCliq] <- PARTITION[cliques[apuCliq]]
  sepPartition[apuSep] <- PARTITION[separators[apuSep]]

  cliqcounts <- zeros(ncliq, npops)
  for (i in 1:npops) {
    cliqcounts[, i] <- rowSums(cliqPartition == i)
  }

  sepcounts <- zeros(nsep, npops)
  for (i in 1:npops) {
    sepcounts[, i] <- rowSums(sepPartition == i)
  }

  return(list(cliqcounts = cliqcounts, sepcounts = sepcounts))
}
