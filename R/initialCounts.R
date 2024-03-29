initialCounts <- function(partition, data, npops, rows, noalle, adjprior) {
  nloci <- size(data, 2)
  ninds <- size(rows, 1)

  koot <- rows[1] - rows[2] + 1
  maxSize <- base::max(koot)

  counts <- zeros(base::max(noalle), nloci, npops)
  sumcounts <- zeros(npops, nloci)
  for (i in 1:npops) {
    for (j in 1:nloci) {
      havainnotLokuksessa <- matlab2r::find(partition == i & data[, j] >= 0)
      sumcounts[i, j] <- length(havainnotLokuksessa)
      for (k in 1:noalle[j]) {
        alleleCode <- k
        N_ijk <- length(
          matlab2r::find(data[havainnotLokuksessa, j] == alleleCode)
        )
        counts[k, j, i] <- N_ijk
      }
    }
  }
  logml <- laskeLoggis(counts, sumcounts, adjprior)
  return(list(sumcounts = sumcounts, counts = counts, logml = logml))
}
