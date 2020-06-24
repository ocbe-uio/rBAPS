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
	merkki <- line[n]
	while (merkki != ',') {
		n <- n + 1
		merkki <- line[n]
	}

	while (!any(merkki == '0123456789')) {
		n <- n + 1
		merkki <- line[n]
	}
	numeroja <- 0
	while (any(merkki == '0123456789')) {
		numeroja <- numeroja + 1
		n <- n + 1
		merkki <- line[n]
	}

	df <- numeroja / 2
	return(df)
}