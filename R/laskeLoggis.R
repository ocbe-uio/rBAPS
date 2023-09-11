laskeLoggis <- function(counts, sumcounts, adjprior) {
  npops <- size(counts, 3)
  replicated_adjprior <- array(adjprior, c(nrow(adjprior), ncol(adjprior), npops))
  sum1 <- sum(sum(sum(lgamma(counts + replicated_adjprior))))
  sum3 <- sum(sum(lgamma(adjprior))) - sum(sum(lgamma(1 + sumcounts)))
  logml2 <- sum1 - npops * sum3
  loggis <- logml2
  return(loggis)
}
