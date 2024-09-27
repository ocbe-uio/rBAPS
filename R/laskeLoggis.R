laskeLoggis <- function(counts, sumcounts, adjprior) {
  npops <- size(counts, 3)
  replicated_adjprior <- array(adjprior, c(nrow(adjprior), ncol(adjprior), npops))
  sum1 <- sum(sum(sum(lgamma(counts + replicated_adjprior))))
  sum2 <- npops * sum(sum(lgamma(adjprior))) + sum(sum(lgamma(1 + sumcounts)))
  sum1 - sum2
}
