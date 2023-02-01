getPopDistancesByKL <- function(adjprior) {
  # Laskee populaatioille etהisyydet
  # kהyttהen KL-divergenssi?
  COUNTS <- COUNTS[seq_len(nrow(adjprior)), seq_len(ncol(adjprior)), ]
  maxnoalle <- size(COUNTS, 1)
  nloci <- size(COUNTS, 2)
  npops <- size(COUNTS, 3)
  distances <- zeros(choose(npops, 2), 1)

  d <- zeros(maxnoalle, nloci, npops)
  prior <- adjprior
  prior[find(prior == 1)] <- 0

  # Lokukset, joissa oli havaittu vain yht?alleelia.
  nollia <- find(all(prior == 0))

  prior[1, nollia] <- 1
  for (pop1 in 1:npops) {
    d[, , pop1] <- (squeeze(COUNTS[, , pop1]) + prior) / repmat(
      sum(squeeze(COUNTS[, , pop1]) + prior), c(maxnoalle, ncol(prior))
    )
  }
  pointer <- 1
  for (pop1 in 1:(npops - 1)) {
    for (pop2 in (pop1 + 1):npops) {
      dist1 <- d[, , pop1]
      dist2 <- d[, , pop2]
      div12 <- sum(sum(dist1 * base::log2((dist1 + 10^-10) / (dist2 + 10^-10)))) /
        nloci
      div21 <- sum(sum(dist2 * base::log2((dist2 + 10^-10) / (dist1 + 10^-10)))) /
        nloci
      div <- (div12 + div21) / 2
      distances[pointer] <- div
      pointer <- pointer + 1
    }
  }
  Z <- linkage(t(distances))
  return(list(Z = Z, distances = distances))
}
