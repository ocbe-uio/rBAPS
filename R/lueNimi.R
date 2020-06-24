#' @title Read the Name
#' @description Reads the line name
#' @param line line
#' @return nimi
#' @export
lueNimi <- function(line) {
	# Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
	n <- 1
	merkki <- line[n]
	nimi <- ''
	while (merkki != ',') {
		nimi <- c(nimi, merkki)
		n <- n + 1
		merkki <- line[n]
	}
	return(nimi)
}