#' @title Handle Pop data
#' @description The last column of the original data tells you which
#' <missing translation> that line is from. The function changes the allele
#' codes so that the codes for one locus have values between 1 and noalle[j].
#' Before this change, an allele whose code is zero is changed.
#' @param raw_data raw data
#' @export
handlePopData <- function(raw_data) {
  # Alkuperäisen datan viimeinen sarake kertoo, milt?yksilölt?
  # kyseinen rivi on peräisin. Funktio muuttaa alleelikoodit
  # siten, ett?yhden lokuksen j koodit saavat arvoja
  # välill?1, ,noalle(j). Ennen tät?muutosta alleeli, jonka
  # koodi on nolla muutetaan.

  data <- raw_data
  nloci <- size(raw_data, 2) - 1

  dataApu <- data[, 1:nloci]
  nollat <- find(dataApu == 0)
  if (length(nollat) > 0) {
     isoinAlleeli <- max(max(dataApu)$maxs)$maxs
     dataApu[nollat] <- isoinAlleeli + 1
     data[, 1:nloci] <- dataApu
  }

  noalle <- zeros(1, nloci)
  alleelitLokuksessa <- cell(nloci, 1, expandable = TRUE)
  for (i in 1:nloci) {
      alleelitLokuksessaI <- sort(unique(data[, i]))
      alleelitLokuksessa[[i]]  <- alleelitLokuksessaI[find(alleelitLokuksessaI >= 0)]
      noalle[i] <- length(alleelitLokuksessa[[i]])
  }
  alleleCodes <- zeros(unique(max(noalle)$maxs), nloci)
  for (i in 1:nloci) {
      alleelitLokuksessaI <- alleelitLokuksessa[[i]]
      puuttuvia <- unique(max(noalle)$maxs) - length(alleelitLokuksessaI)
      alleleCodes[, i] = c(alleelitLokuksessaI, zeros(puuttuvia, 1))
  }

  for (loc in 1:nloci) {
      for (all in 1:noalle[loc]) {
          data[find(data[ , loc] == alleleCodes[all, loc]), loc] <- all
      }
  }

  nind <- max(data[, ncol(data)])$maxs
  rows <- zeros(nind, 2)
  for (i in 1:nind) {
      rivit <- t(find(data[, ncol(data)] == i))
      rows[i, 1] <- min(rivit)$mins
      rows[i, 2] <- max(rivit)$maxs
  }
  newData <- data

  adjprior <- zeros(unique(max(noalle)$maxs), nloci)
  priorTerm <- 0
  for (j in 1:nloci) {
      adjprior[, j]  <- c(repmat(1 / noalle[j], c(noalle[j], 1)), ones(unique(max(noalle)$maxs) - noalle[j], 1))
      priorTerm <- priorTerm + noalle[j] * log(gamma(1 / noalle[j]))
  }
	return(
        list(
            newData = newData,
            rows = rows,
            alleleCodes = alleleCodes,
            noalle = noalle,
            adjprior = adjprior,
            priorTerm = priorTerm
        )
    )
}
