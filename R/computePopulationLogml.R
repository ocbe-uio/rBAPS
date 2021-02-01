computePopulationLogml <- function(pops, adjprior, priorTerm = 0) {
	# Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset

	# ======================================================== #
	# Limiting COUNTS size                                     #
	# ======================================================== #
	COUNTS <- COUNTS[seq_len(nrow(adjprior)), seq_len(ncol(adjprior)), pops, drop=FALSE]

	x <- size(COUNTS, 1)
	y <- size(COUNTS, 2)
	z <- length(pops)

	# ======================================================== #
	# Computation                                              #
	# ======================================================== #
	term1 <- squeeze(
		# FIXME: assumes COUNTS has 3 dims. Where does this come from?
		sum(
			sum(
				reshape(
					lgamma(
						repmat(adjprior, c(1, 1, length(pops))) +
							COUNTS[seq_len(nrow(adjprior)), seq_len(ncol(adjprior)), pops, drop=FALSE]
					),
					c(x, y, z)
				),
				1
			),
			2
		)
	)
	popLogml <- term1 - sum(lgamma(1 + SUMCOUNTS[pops, ]), 2) - priorTerm
	return(popLogml)
}