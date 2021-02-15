updateGlobalVariables <- function(ind, i2, diffInCounts, adjprior, priorTerm) {
	# % Suorittaa globaalien muuttujien muutokset, kun yksil� ind
	# % on siirret��n koriin i2.
	i1 <- PARTITION[ind]
	PARTITION[ind] <- i2

	COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
	COUNTS[, , i2] <- COUNTS[, , i2] + diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - colSums(diffInCounts)
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + colSums(diffInCounts)

	POP_LOGML[c(i1, i2)] <- computePopulationLogml(
		c(i1, i2), adjprior, priorTerm
	)

	LOGDIFF[, c(i1, i2)] <- -Inf
	inx <- c(find(PARTITION == i1), find(PARTITION==i2))
	LOGDIFF[inx, ] <- -Inf
}

updateGlobalVariables2 <- function(i1, i2, diffInCounts, adjprior, priorTerm) {
	# % Suorittaa globaalien muuttujien muutokset, kun kaikki
	# % korissa i1 olevat yksil�t siirret��n koriin i2.

	inds <- find(PARTITION == i1)
	PARTITION[inds] <- i2

	COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
	COUNTS[, , i2] <- COUNTS[, , i2] + diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - colSums(diffInCounts)
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + colSums(diffInCounts)

	POP_LOGML[i1] <- 0
	POP_LOGML[i2] <- computePopulationLogml(i2, adjprior, priorTerm)

	LOGDIFF[, c(i1, i2)] <- -Inf
	inx <- c(find(PARTITION == i1), find(PARTITION == i2))
	LOGDIFF[inx, ] <- -Inf
}

updateGlobalVariables3 <- function(
	muuttuvat, diffInCounts, adjprior, priorTerm, i2
) {
	# % Suorittaa globaalien muuttujien p�ivitykset, kun yksil�t 'muuttuvat'
	# % siirret��n koriin i2. Ennen siirtoa yksil�iden on kuuluttava samaan
	# % koriin.

	i1 <- PARTITION[muuttuvat(1)]
	PARTITION[muuttuvat] <- i2

	COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
	COUNTS[, , i2] <- COUNTS[, , i2] + diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - colSums(diffInCounts)
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + colSums(diffInCounts)

	POP_LOGML[c(i1, i2)] <- computePopulationLogml(
		c(i1, i2), adjprior, priorTerm
	)

	LOGDIFF[, c(i1, i2)] <- -Inf
	inx <- c(find(PARTITION == i1), find(PARTITION == i2))
	LOGDIFF[inx, ] <- -Inf
}
