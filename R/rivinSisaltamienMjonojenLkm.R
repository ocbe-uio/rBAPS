#' @title Number of M queues
#' @param line line number
#' @return count
#' @export
rivinSisaltamienMjonojenLkm <- function(line) {
	# Palauttaa line:n sis�lt�mien mjonojen lukum��r�n.
	# Mjonojen v�liss?t�ytyy olla v�lily�nti.
	count <- 0
	pit <- length(line)
	tila <- 0 # 0, jos odotetaan v�lily�ntej? 1 jos odotetaan muita merkkej?
	for (i in seq_len(pit)) {
		merkki <- line[i]
		if (merkki == " " & tila == 0) {
			# Ei tehd?mit��n.
		} else if (merkki == " " & tila == 1) {
			tila <- 0
		} else if (merkki != " " & tila == 0) {
			tila <- 1
			count <- count + 1
		} else if (merkki != " " & tila == 1) {
			# %Ei tehd?mit��n
		}
	}
	return(count)
}