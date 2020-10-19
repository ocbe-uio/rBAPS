addToSummary <- funciton(logml, partitionSummary, worstIndex) {
	# Tiedet��n, ett� annettu logml on isompi kuin huonoin arvo
	# partitionSummary taulukossa. Jos partitionSummary:ss� ei viel� ole
	# annettua logml arvoa, niin lis�t��n worstIndex:in kohtaan uusi logml ja
	# nykyist� partitiota vastaava nclusters:in arvo. Muutoin ei tehd� mit��n.

	apu <- find(abs(partitionSummary[, 2] - logml) < 1e-5)
	if (isempty(apu)) {
	    # Nyt l�ydetty partitio ei ole viel� kirjattuna summaryyn.
	    npops <- length(unique(PARTITION))
	    partitionSummary[worstIndex, 1] <- npops
	    partitionSummary[worstIndex, 2] <- logml
	    added <- 1
	} else {
	    added <- 0
	}
	return(list(partitionSummary = partitionSummary, added = added))
}