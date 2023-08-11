#' @title Handle Data
#' @param raw_data Raw data in Genepop or BAPS format
#' @param format data format
#' @details The last column of the original data tells you from which
#' individual that line is from. The function first examines how many line
#' maximum is from one individual giving know if it is haploid, diploid, etc.
#' After this function. Add blank lines for individuals with fewer rows as
#' maximum. If the code of an allele is = 0, the function changes that allele
#' code to the smallest code that is larger than any code in use. After this,
#' the function changes the allele codes so that one locus j
#' codes get values between? 1, ..., noalle(j).
handleData <- function(raw_data, format = "Genepop") {
  # Alkuper?isen datan viimeinen sarake kertoo, milt?yksil?lt?
  # kyseinen rivi on per?isin. Funktio tutkii ensin, ett?montako
  # rivi?maksimissaan on per?isin yhdelt?yksil?lt? jolloin saadaan
  # tiet?? onko kyseess?haploidi, diploidi jne... T?m?n j?lkeen funktio
  # lis?? tyhji?rivej?niille yksil?ille, joilta on per?isin v?hemm?n
  # rivej?kuin maksimim??r?
  #   Mik?li jonkin alleelin koodi on =0, funktio muuttaa t?m?n alleelin
  # koodi pienimm?ksi koodiksi, joka isompi kuin mik??n k?yt?ss?oleva koodi.
  # T?m?n j?lkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
  # koodit saavat arvoja v?lill?1,...,noalle(j).
  data <- as.matrix(raw_data)
  if (format %in% c("genepop", "baps")) {
    nloci <- size(raw_data, 2) - 1
  } else {
    nloci <- size(raw_data, 2)
  }

  dataApu <- data[, 1:nloci]
  nollat <- matlab2r::find(dataApu == 0)
  if (!isempty(nollat)) {
    isoinAlleeli <- base::max(base::max(dataApu))
    dataApu[nollat] <- isoinAlleeli + 1
    data[, 1:nloci] <- dataApu
  }

  noalle <- zeros(1, nloci)
  alleelitLokuksessa <- cell(nloci, 1, expandable = TRUE)
  for (i in 1:nloci) {
    alleelitLokuksessaI <- unique(data[, i])
    alleelitLokuksessa[[i]] <- sort(alleelitLokuksessaI[
      matlab2r::find(alleelitLokuksessaI >= 0)
    ])
    noalle[i] <- length(alleelitLokuksessa[[i]])
  }
  alleleCodes <- zeros(base::max(noalle), nloci)
  for (i in 1:nloci) {
    alleelitLokuksessaI <- alleelitLokuksessa[[i]]
    puuttuvia <- base::max(noalle) - length(alleelitLokuksessaI)
    alleleCodes[, i] <- as.matrix(c(alleelitLokuksessaI, zeros(puuttuvia, 1)))
  }

  for (loc in seq_len(nloci)) {
    for (all in seq_len(noalle[loc])) {
      data[matlab2r::find(data[, loc] == alleleCodes[all, loc]), loc] <- all
    }
  }

  nind <- as.integer(base::max(data[, ncol(data)]))
  nrows <- size(data, 1)
  ncols <- size(data, 2)
  rowsFromInd <- zeros(nind, 1)
  for (i in 1:nind) {
    rowsFromInd[i] <- length(matlab2r::find(data[, ncol(data)] == i))
  }
  maxRowsFromInd <- base::max(rowsFromInd)
  a <- -999
  emptyRow <- repmat(a, c(1, ncols))
  lessThanMax <- matlab2r::find(rowsFromInd < maxRowsFromInd)
  missingRows <- max(maxRowsFromInd * nind - nrows, 0L)
  data <- rbind(data, zeros(missingRows, ncols))
  pointer <- 1
  for (ind in t(lessThanMax)) { # K?y l?pi ne yksil?t, joilta puuttuu rivej?
    miss <- maxRowsFromInd - rowsFromInd[ind] # T?lt?yksil?lt?puuttuvien lkm.
  }
  data <- sortrows(data, ncols) # Sorttaa yksil?iden mukaisesti
  newData <- data
  rowsFromInd <- maxRowsFromInd

  adjprior <- zeros(base::max(noalle), nloci)
  priorTerm <- 0
  for (j in 1:nloci) {
    adjprior[, j] <- as.matrix(c(
      repmat(1 / noalle[j], c(noalle[j], 1)),
      ones(base::max(noalle) - noalle[j], 1)
    ))
    priorTerm <- priorTerm + noalle[j] * lgamma(1 / noalle[j])
  }
  out <- list(
    newData     = newData,
    rowsFromInd = rowsFromInd,
    alleleCodes = alleleCodes,
    noalle      = noalle,
    adjprior    = adjprior,
    priorTerm   = priorTerm
  )
  return(out)
}
