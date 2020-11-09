initialPopCounts <- function(data, npops, rows, noalle, adjprior) {
	nloci <- size(data, 2)
	counts <- zeros(max(noalle), nloci, npops)
	sumcounts <- zeros(npops, nloci)

	for (i in 1:npops) {
		for (j in 1:nloci) {
			i_rivit <- rows(i, 1):rows(i, 2)
			havainnotLokuksessa <- find(data[i_rivit, j] >= 0)
			sumcounts(i, j) <- length(havainnotLokuksessa)
			for (k in 1:noalle[j]) {
				alleleCode <- k
				N_ijk <- length(find(data[i_rivit, j] == alleleCode))
				counts(k, j, i) <- N_ijk
			}
		}
	}
	logml <- laskeLoggis(counts, sumcounts, adjprior)
	return(sumcounts = sumcounts, counts = counts, logml = logml)
}
