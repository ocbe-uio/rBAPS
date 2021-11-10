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

  nloci <- size(data, 2) # added 1 from original code
  if (size(data, 1) < (2 * ind)) {
    data <- rbind(data, zeros(100, nloci)) # subtracted 1 from original code
  }

  k <- 1
  merkki <- substring(line, k, k)
  while (merkki != ",") {
    k <- k + 1
    merkki <- substring(line, k, k)
  }
  line <- substring(line, k + 1)
  # clear k; clear merkki;

  if (grepl(" ", line)) {
    alleeliTaulu <- as.numeric(strsplit(line, split = " ")[[1]])
  } else if (grepl("\t", line)) {
    alleeliTaulu <- as.numeric(strsplit(line, split = "\t")[[1]])
  }

  if (length(alleeliTaulu) != nloci) {
    stop("Incorrect data format.")
  }

  for (j in seq_len(nloci)) {
    ekaAlleeli <- floor(alleeliTaulu[j] / divider)
    if (is.na(ekaAlleeli) | ekaAlleeli == 0) ekaAlleeli <- -999
    tokaAlleeli <- alleeliTaulu[j] %% divider
    if (is.na(tokaAlleeli) | tokaAlleeli == 0) tokaAlleeli <- -999

    data[2 * ind - 1, j] <- ekaAlleeli
    data[2 * ind, j] <- tokaAlleeli
  }

  data[2 * ind - 1, ncol(data)] <- ind
  data[2 * ind, ncol(data)] <- ind
  return(data)
}
