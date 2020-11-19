findEmptyPop <- function(npops) {
	# % Palauttaa ensimm�isen tyhj�n populaation indeksin. Jos tyhji�
	# % populaatioita ei ole, palauttaa -1:n.
	pops <- t(unique(PARTITION))
	if (length(pops) == npops) {
		emptyPop <- -1
	} else {
		popDiff <- diff(c(0, pops, npops + 1))
		emptyPop <- min(find(popDiff > 1))
	}
	return(list(emptyPop = emptyPop, pops = pops))
}
