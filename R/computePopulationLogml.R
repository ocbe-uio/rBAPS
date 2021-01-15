computePopulationLogml <- function(pops, adjprior, priorTerm) {
	# Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset

	x <- size(COUNTS, 1)
	y <- size(COUNTS, 2)
	z <- length(pops)

	popLogml <- squeeze(
		# FIXME: assumes COUNTS has 3 dims. Where does this come from?
		sum(
			sum(
				reshape(
					lgamma(
						repmat(adjprior, c(1, 1, length(pops))) +
							COUNTS[, , pops]
					),
					c(x, y, z)
				),
				1
			),
			2
		)
	) - sum(lgamma(1 + SUMCOUNTS[pops, ]), 2) - priorTerm
	return(popLogml)
}