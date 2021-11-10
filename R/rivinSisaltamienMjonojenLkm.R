#' @title Number of M queues
#' @param line line number
#' @return count
#' @description Returns the number of queues contained in the line. There must be a space between the queues.
#' @export
rivinSisaltamienMjonojenLkm <- function(line) {
  # Palauttaa line:n sis�lt�mien mjonojen lukum��r�n.
  # Mjonojen v�liss?t�ytyy olla v�lily�nti.
  count <- 0
  pit <- nchar(line)
  tila <- 0 # 0, jos odotetaan v�lily�ntej? 1 jos odotetaan muita merkkej?
  for (i in seq_len(pit)) {
    merkki <- substring(line, i, i)
    if (isspace(merkki) & tila == 0) {
      # Ei tehd?mit��n.
    } else if (isspace(merkki) & tila == 1) {
      tila <- 0
    } else if (!isspace(merkki) & tila == 0) {
      tila <- 1
      count <- count + 1
    } else if (!isspace(merkki) & tila == 1) {
      # %Ei tehd?mit��n
    }
  }
  return(count)
}
