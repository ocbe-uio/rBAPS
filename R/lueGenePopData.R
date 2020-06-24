#' @title Read GenePop Data
#' @description Reads GenePop-formatted data
#' @param tiedostonNimi Name of the file
#' @return list containing data and popnames
#' @export
lueGenePopData <- function (tiedostonNimi) {

	fid   <- load(tiedostonNimi)
	line1 <- readLines(fid)[1] # ensimmäinen rivi
	line2 <- readLines(fid)[2] # toinen rivi
	count <- rivinSisaltamienMjonojenLkm(line)

	line <- readLines(fid)[3]
	lokusRiveja <- 1
	while (testaaPop(line) == 0) {
		lokusRiveja <- lokusRiveja + 1 # locus row
		line <- readLines(fid)[3 + lokusRiveja]
	}

	if (lokusRiveja > 1) {
		nloci <- lokusRiveja
	} else {
		nloci <- count
	}

	popnames <- cell(10, 2)
	data <- zeros(100, nloci + 1)
	nimienLkm <- 0
	ninds <- 0
	poimiNimi <- 1
	digitFormat <- -1
	while (line != -1) {
		line <- readLines(fid)[lokusRiveja + 1]
		lokusRiveja <- lokusRiveja + 1

		if (poimiNimi == 1) {
			# Edellinen rivi oli 'pop'
			nimienLkm <- nimienLkm + 1
			ninds <- ninds + 1
			if (nimienLkm > size(popnames, 1)) {
				popnames <- c(popnames, cell(10, 2))
			}
			nimi <- lueNimi(line)
			if (digitFormat == -1) {
				digitFormat <- selvitaDigitFormat(line)
				divider <- 10 ^ digitFormat
			}
			popnames[nimienLkm, 1] <- nimi   #N�in se on greedyMix:iss�kin?!?
			popnames[nimienLkm, 2] <- ninds
			poimiNimi <- 0

			data <- addAlleles(data, ninds, line, divider)

		 } else if (testaaPop(line)) {
			poimiNimi <- 1

		 } else if (line != -1) {
			ninds <- ninds + 1
			data <- addAlleles(data, ninds, line, divider)
		}
	}

	data <- data[1:(ninds * 2),]
	popnames <- popnames[seq_len(nimienLkm),]
	return(list(data = data, popnames = popnames))
}