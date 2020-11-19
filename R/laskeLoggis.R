laskeLoggis <- function(counts, sumcounts, adjprior) {
	npops <- size(counts, 3)

	sum1 <- sum(sum(sum(lgamma(counts + repmat(adjprior, c(1, 1, npops))))))
	sum3 <- sum(sum(lgamma(adjprior))) - sum(sum(lgamma(1 + sumcounts)))
	logml2 <- sum1 - npops * sum3
	loggis <- logml2
	return(loggis)
}