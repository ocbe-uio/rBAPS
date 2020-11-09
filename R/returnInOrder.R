returnInOrder <- function(inds, pop, globalRows, data, adjprior, priorTerm) {
	# % Palauttaa yksil�t j�rjestyksess� siten, ett� ensimm�isen� on
	# % se, jonka poistaminen populaatiosta pop nostaisi logml:n
	# % arvoa eniten.

	ninds <- length(inds)
	apuTaulu <- [inds, zeros(ninds,1)];

	for (i in 1:ninds) {
		ind <- inds[i]
		rows <- globalRows[i, 1]:globalRows[i, 2]
		diffInCounts <- computeDiffInCounts(
			rows, size[COUNTS, 1], size[COUNTS, 2], data
		)
		diffInSumCounts <- sum(diffInCounts)

		COUNTS[ , ,pop] <- COUNTS[ , ,pop] - diffInCounts
		SUMCOUNTS[pop, ] <- SUMCOUNTS[pop, ] - diffInSumCounts
		apuTaulu[i, 2] <- computePopulationLogml(pop, adjprior, priorTerm)
		COUNTS[ , ,pop] <- COUNTS[ , ,pop] + diffInCounts
		SUMCOUNTS[pop, ] <- SUMCOUNTS[pop, ] + diffInSumCounts
	}
	apuTaulu <- sortrows(apuTaulu, 2)
	inds <- apuTaulu[ninds:1, 1]
	return(inds)
}