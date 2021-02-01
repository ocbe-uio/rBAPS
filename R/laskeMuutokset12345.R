#' @title Calculate changes (?)
#' @description Palauttaa npops*npops taulun, jonka alkio (i,j) kertoo, mik?on
#' muutos logml:ss? mikäli populaatiosta i siirretään osuuden verran
#' todennäköisyysmassaa populaatioon j. Mikäli populaatiossa i ei ole mitään
#' siirrettävää, on vastaavassa kohdassa rivi nollia.
#' @param osuus Percentages?
#' @param omaFreqs own Freqs?
#' @param osuusTaulu Percentage table?
#' @param logml log maximum likelihood
#' @export
laskeMuutokset4 <- function (osuus, osuusTaulu, omaFreqs, logml) {
	if (is.null(dim(COUNTS))) {
		npops <- 1
	} else {
		npops <- ifelse(is.na(dim(COUNTS)[3]), 1, dim(COUNTS)[3])
	}
	notEmpty <- which(osuusTaulu > 0.005)
	muutokset <- zeros(npops)
	empties <- !notEmpty

	for (i1 in notEmpty) {
		osuusTaulu[i1] <- osuusTaulu[i1] - osuus
		for (i2 in c(colon(1, i1 - 1), colon(i1 + 1, npops))) {
			osuusTaulu[i2] <- osuusTaulu[i2] + osuus
			loggis <- computeIndLogml(omaFreqs, osuusTaulu)

			# Work around Matlab OOB bug
			if (i1 > nrow(muutokset)) {
				muutokset <- rbind(muutokset, muutokset * 0)
			}
			if (i2 > ncol(muutokset)) {
				muutokset <- cbind(muutokset, muutokset * 0)
			}

			muutokset[i1, i2] <- loggis - logml
			osuusTaulu[i2] <- osuusTaulu[i2] - osuus
		}
		osuusTaulu[i1] <- osuusTaulu[i1] + osuus
	}
	return (muutokset)
}

# Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik� olisi
# muutos logml:ss�, mik�li yksil� ind siirret��n koriin i.
# diffInCounts on poistettava COUNTS:in siivusta i1 ja lis�tt�v�
# COUNTS:in siivuun i2, mik�li muutos toteutetaan.
#
# Lis�ys 25.9.2007:
# Otettu k�ytt��n globaali muuttuja LOGDIFF, johon on tallennettu muutokset
# logml:ss� siirrett�ess� yksil�it� toisiin populaatioihin.
laskeMuutokset <- function(ind, globalRows, data, adjprior, priorTerm) {
	npops <- size(COUNTS, 3)
	muutokset <- LOGDIFF[ind, ]

	i1 <- PARTITION[ind]
	i1_logml <- POP_LOGML[i1]
	muutokset[i1] <- 0

	rows <- globalRows[ind, 1]:globalRows[ind, 2]
	diffInCounts <- computeDiffInCounts(
		rows, size(COUNTS, 1), size(COUNTS, 2), data
	)
	diffInSumCounts <- sum(diffInCounts)

	COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
	new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
	COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

	# TODO: check i2 calculation against MATLAB (where does this code come from?)
	i2 <- find(muutokset == -Inf) # Etsit��n populaatiot jotka muuttuneet viime kerran j�lkeen. (Searching for populations that have changed since the last time)
	i2 <- setdiff(i2, i1)
	i2_logml <- POP_LOGML[i2]

	ni2 <- length(i2)

	# FIXME: i2 is empty
	browser() # TEMP
	COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, ni2))
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(ni2, 1))
	new_i2_logml <- computePopulationLogml(i2, adjprior, priorTerm)
	COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, ni2))
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(ni2, 1))

	muutokset[i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
	LOGDIFF[ind, ] = muutokset
	return(list(muutokset = muutokset, diffInCounts = diffInCounts))
}

laskeMuutokset2 <- function(i1, globalRows, data, adjprior, priorTerm) {
	# % Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik� olisi
	# % muutos logml:ss�, mik�li korin i1 kaikki yksil�t siirret��n
	# % koriin i.

	npops <- size(COUNTS, 3)
	muutokset <- zeros(npops, 1)

	i1_logml <- POP_LOGML[i1]

	inds <- find(PARTITION == i1)
	ninds <- length(inds)

	if (ninds == 0) {
		diffInCounts <- zeros(size(COUNTS, 1), size(COUNTS, 2))
		return()
	}

	rows = list()
	for (i in 1:ninds) {
		ind <- inds(i)
		lisa <- globalRows(ind, 1):globalRows(ind, 2)
		rows <- c(rows, t(lisa))
	}

	diffInCounts <- computeDiffInCounts(
		t(rows), size(COUNTS, 1), size(COUNTS, 2), data
	)
	diffInSumCounts <- sum(diffInCounts)

	COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
	new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
	COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
	SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

	i2 <- c(1:i1-1, i1+1:npops)
	i2_logml <- POP_LOGML[i2]

	COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, npops - 1))
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(npops - 1, 1))
	new_i2_logml <- computePopulationLogml(i2, adjprior, priorTerm)
	COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, npops - 1))
	SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(npops - 1, 1))

	muutokset[i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
	return(list(muutokset = muutokset, diffInCounts = diffInCounts))
}


laskeMuutokset3 <- function(T2, inds2, globalRows, data, adjprior, priorTerm, i1) {
	# Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
	# kertoo, mik� olisi muutos logml:ss�, jos populaation i1 osapopulaatio
	# inds2(find(T2==i)) siirret��n koriin j.

	npops <- size(COUNTS, 3)
	npops2 <- length(unique(T2))
	muutokset <- zeros(npops2, npops)

	i1_logml = POP_LOGML[i1]
	for (pop2 in 1:npops2) {
		inds <- inds2[find(T2==pop2)]
		ninds <- length(inds);
		if (ninds > 0) {
			rows <- list()
			for (i in 1:ninds) {
				ind <- inds[i]
				lisa <- globalRows[ind, 1]:globalRows[ind, 2]
				rows <- c(rows, t(lisa))
			}
			diffInCounts <- computeDiffInCounts(
				t(rows), size(COUNTS, 1), size(COUNTS, 2), data
			)
			diffInSumCounts <- sum(diffInCounts)

			COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
			SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
			new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
			COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
			SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

			i2 <- c(1:i1-1, i1+1:npops)
			i2_logml <- t(POP_LOGML[i2])

			COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, npops - 1))
			SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(npops - 1, 1))
			new_i2_logml <- t(computePopulationLogml(i2, adjprior, priorTerm))
			COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, npops - 1))
			SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(npops - 1, 1))

			muutokset[pop2, i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
		}
	}
	return(muutokset)
}

laskeMuutokset5 <- function(inds, globalRows, data, adjprior, priorTerm, i1, i2) {
	# Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik� olisi
	# muutos logml:ss�, mik�li yksil� i vaihtaisi koria i1:n ja i2:n v�lill�.

	ninds <- length(inds)
	muutokset <- zeros(ninds, 1)

	i1_logml <- POP_LOGML[i1]
	i2_logml <- POP_LOGML[i2]

	for (i in 1:ninds) {
		ind <- inds[i]
		if (PARTITION[ind] == i1) {
			pop1 <- i1  #mist�
			pop2 <- i2  #mihin
		} else {
			pop1 <- i2
			pop2 <- i1
		}
		rows <- globalRows[ind, 1]:globalRows[ind, 2]
		diffInCounts <- computeDiffInCounts(
			rows, size(COUNTS, 1), size(COUNTS, 2), data
		)
		diffInSumCounts <- sum(diffInCounts)



		COUNTS[, , pop1] <- COUNTS[, , pop1] - diffInCounts
		SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] - diffInSumCounts
		COUNTS[, , pop2] <- COUNTS[, , pop2] + diffInCounts
		SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] + diffInSumCounts

		new_logmls <- computePopulationLogml(c(i1, i2), adjprior, priorTerm)
		muutokset[i] <- sum(new_logmls)

		COUNTS[, , pop1] <- COUNTS[, , pop1] + diffInCounts
		SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] + diffInSumCounts
		COUNTS[, , pop2] <- COUNTS[, , pop2] - diffInCounts
		SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] - diffInSumCounts
	}

	muutokset <- muutokset - i1_logml - i2_logml
	return(muutokset)
}