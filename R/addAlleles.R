#' @title Add Alleles
#' @param data data
#' @param ind ind
#' @param line line
#' @param divider divider
#' @return data (after alleles were added)
#' @export
addAlleles <- function(data, ind, line, divider) {
	# Lisaa BAPS-formaatissa olevaan datataulukkoon
	# yksil�� ind vastaavat rivit. Yksil�n alleelit
	# luetaan genepop-formaatissa olevasta rivist?
	# line. Jos data on 3 digit formaatissa on divider=1000.
	# Jos data on 2 digit formaatissa on divider=100.

	nloci <- size(data, 2) - 1
	if (size(data, 1) < (2 * ind)) {
		data <- c(data, zeros(100, nloci + 1))
	}

	k <- 1
	merkki <- line[k]
	while (merkki != ',') {
	   k <- k + 1
	   merkki <- line[k]
	}
	line <- line[k + 1:length(line)]
	# clear k; clear merkki;

	alleeliTaulu <- as.numeric(strsplit(line, split = " ")[[1]])


	if (length(alleeliTaulu) != nloci) {
		stop('Incorrect data format.')
	}

	for (j in seq_len(nloci)) {
		ekaAlleeli <- floor(alleeliTaulu[j] / divider)
		if (ekaAlleeli == 0) ekaAlleeli <- -999
		tokaAlleeli <- alleeliTaulu[j] %% divider
		if (tokaAlleeli == 0) tokaAlleeli <- -999

		data[2 * ind - 1, j] <- ekaAlleeli
		data[2 * ind, j] <- tokaAlleeli
	}

	data[2 * ind - 1,end] <- ind
	data[2 * ind, end] <- ind
	return(data)
}