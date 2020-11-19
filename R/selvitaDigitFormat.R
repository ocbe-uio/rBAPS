#' @title Find out the Digit Format
#' @param line the first line after the "pop" word from data in Genepop format. #  @note Function clarified based on the line format whether the alleles of the data are given using 2 or 3 numbers.
#' @return df
#' @export
selvitaDigitFormat <- function(line) {
	# line on ensimm�inen pop-sanan j�lkeinen rivi
	# Genepop-formaatissa olevasta datasta. funktio selvitt��
	# rivin muodon perusteella, ovatko datan alleelit annettu
	# 2 vai 3 numeron avulla.
	n <- 1
	merkki <- substring(line, n, n)
	while (merkki != ',') {
		n <- n + 1
		merkki <- substring(line, n, n)
	}

	while (!any(merkki %in% as.character(0:9))) {
		n <- n + 1
		merkki <- substring(line, n, n)
	}
	numeroja <- 0
	while (any(merkki %in% as.character(0:9))) {
		numeroja <- numeroja + 1
		n <- n + 1
		merkki <- substring(line, n, n)
	}

	df <- numeroja / 2
	return(df)
}