lueGenePopDataPop <- function(tiedostonNimi) {
  # Data annetaan muodossa, jossa viimeinen sarake kertoo ryhmän.
  # popnames on kuten ennenkin.

  fid <- fopen(tiedostonNimi)
  line <- fgetl(fid)  # ensimmäinen rivi
  line <- fgetl(fid)  # toinen rivi
  count <- rivinSisaltamienMjonojenLkm(line)

  line <- fgetl(fid)
  lokusRiveja <- 1
  while (testaaPop(line) == 0) {
      lokusRiveja <- lokusRiveja + 1
      line <- fgetl(fid)
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
  digitFormat = -1
  while (line != -1) {
      line <- fgetl(fid)

      if (poimiNimi == 1) {
          # Edellinen rivi oli 'pop'
          nimienLkm <- nimienLkm + 1
          ninds <- ninds + 1
          if (nimienLkm > size(popnames, 1)) {
              popnames <- c(popnames, cell(10, 2))
          }
          nimi <- lueNimi(line)
          if (digitFormat == -1) {
              digitFormat <- selvitaDigitFormat(line)
              divider <- 10 ^ digitFormat
          }
          popnames[nimienLkm, 1] = nimi   # Näin se on greedyMix:issäkin?!?
          popnames[nimienLkm, 2] = ninds
          poimiNimi <- 0

          data <- addAlleles(data, ninds, line, divider)

      } else if (testaaPop(line)) {
          poimiNimi <- 1

      } else if (line != -1) {
          ninds <- ninds + 1
          data <- addAlleles(data, ninds, line, divider)
      }
  }

  data <- data[1:ninds * 2, ]
  popnames <- popnames[1:nimienLkm, ]
  npops <- size(popnames, 1)
  ind <- 1
  for (pop in 1:npops) {
      if (pop < npops) {
          while (ind < popnames[pop + 1, 2]) {
              data[c(ind * 2 - 1, ind * 2), ncol(data)] <- pop
              ind <- ind + 1
          }
      } else {
          while (ind <= ninds) {
              data[c(ind * 2 - 1, ind * 2), ncol(data)] <- pop
              ind <- ind + 1
          }
      }
  }
	return(list(data = data, popnames = popnames))
}
