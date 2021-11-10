#' @title Read GenePop Data
#' @description Reads GenePop-formatted data
#' @param tiedostonNimi Name of the file
#' @return list containing data and popnames
#' @export
lueGenePopData <- function(tiedostonNimi) {
  fid <- readLines(tiedostonNimi)
  line <- fid[1] # ensimmäinen rivi
  line <- fid[2] # toinen rivi
  count <- rivinSisaltamienMjonojenLkm(line)

  line <- fid[3]
  lokusRiveja <- 1
  while (testaaPop(line) == 0) {
    lokusRiveja <- lokusRiveja + 1 # locus row
    line <- fid[2 + lokusRiveja]
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
  while (lokusRiveja < length(fid) - 2) {
    lokusRiveja <- lokusRiveja + 1 # Keeps the loop moving along
    line <- fid[lokusRiveja + 2]
    if (poimiNimi == 1) {
      # Edellinen rivi oli 'pop' (previous line was pop)
      nimienLkm <- nimienLkm + 1
      ninds <- ninds + 1
      if (nimienLkm > size(popnames, 1)) {
        popnames <- rbind(popnames, cell(10, 2))
      }
      nimi <- lueNimi(line)
      if (digitFormat == -1) {
        digitFormat <- selvitaDigitFormat(line)
        divider <- 10^digitFormat
      }
      popnames[nimienLkm, 1] <- nimi # N�in se on greedyMix:iss�kin?!?
      popnames[nimienLkm, 2] <- ninds
      poimiNimi <- 0
      data <- addAlleles(data, ninds, line, divider)
    } else if (testaaPop(line)) {
      poimiNimi <- 1
    } else if (!is.na(line)) {
      ninds <- ninds + 1
      data <- addAlleles(data, ninds, line, divider)
    }
  }

  data <- data[1:(ninds * 2), ]
  popnames <- popnames[seq_len(nimienLkm), ]
  return(list(data = data, popnames = popnames))
}
