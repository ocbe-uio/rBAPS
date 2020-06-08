#' @title Tests GenePop data
testaaGenePopData <- function(tiedostonNimi) {
	# kunnossa == 0, jos data ei ole kelvollinen genePop data.
	# Muussa tapauksessa kunnossa == 1.

	kunnossa <- 0
	if (file.exists(paste0(tiedostonNimi, ".rda"))) {
		fid   <- load(tiedostonNimi)
		line1 <- readLines(fid)[1] # ensimmäinen rivi
		line2 <- readLines(fid)[2] # toinen rivi
		line3 <- readLines(fid)[3] # kolmas
	} else {
		fid <- line1 <- line2 <- line3 <- -1
	}

	if (line1 == -1 | line2 == -1 | line3 == -1) {
		stop('Incorrect file format 1168')
	}
	if (testaaPop(line1) == 1 | testaaPop(line2) == 1) { # TODO: translate function
		stop('Incorrect file format 1172')
	}
	if (testaaPop(line3) == 1) {
		# 2 rivi t�ll�in lokusrivi
		nloci <- rivinSisaltamienMjonojenLkm(line2) # TODO: translate function
		line4 <- readLines(fid)[4]
		if (line4 == -1) stop('Incorrect file format 1180')
		if (!grepl(',', line4)) {
			# Rivin nelj?t�ytyy sis�lt�� pilkku.
			stop('Incorrect file format 1185')
		}
		pointer <- 1
		while (line4[pointer] != ',') { # Tiedet��n, ett?pys�htyy
			pointer <- pointer + 1
		}
		line4 <- line4[(pointer + 1):nchar(line4)] # pilkun j�lkeinen osa
		nloci2 <- rivinSisaltamienMjonojenLkm(line4)
		if (nloci2 != nloci) stop('Incorrect file format 1195')
	} else {
		line <- readLines(fid)[4]
		lineNumb <- 4
		while (testaaPop(line) != 1 & line != -1) {
			line <- readLines(fid)[lineNumb]
			lineNumb <- lineNumb + 1
		}
		if (line == -1) stop('Incorrect file format 1206')
		nloci <- lineNumb - 2
		line4 <- readLines(fid)[4] # Eka rivi pop sanan j�lkeen
		if (line4 == -1) stop('Incorrect file format 1212')
		if (!grepl(',', line4)) {
			# Rivin t�ytyy sis�lt�� pilkku.
			stop('Incorrect file format 1217')
		}
		pointer <- 1
		while (line4[pointer] != ',') { # Tiedet��n, ett?pys�htyy
			pointer <- pointer + 1
		}
		line4 <- line4[(pointer + 1):nchar(line4)] # pilkun j�lkeinen osa
		nloci2 <- rivinSisaltamienMjonojenLkm(line4)
		if (nloci2 != nloci) stop('Incorrect file format 1228')
	}
	kunnossa <- 1
	return(kunnossa)
}