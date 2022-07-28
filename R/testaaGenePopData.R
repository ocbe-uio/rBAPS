#' @title Tests GenePop data
#' @param tiedostonNimi Filename
#' @return kunnossa (binary "ok" condition value) == 0 if the data is not valid
#' genePop data. Otherwise, kunnossa == 1.
#' @details GenePop data are textfiles that follow the GenePop format. This
#' function checks if such file is properly formatted as GenePop.
testaaGenePopData <- function(tiedostonNimi) {
  # kunnossa == 0, jos data ei ole kelvollinen genePop data.
  # Muussa tapauksessa kunnossa == 1.

  kunnossa <- 0
  if (file.exists(tiedostonNimi)) {
    fid <- readLines(tiedostonNimi)
    line1 <- fid[1] # ensimmäinen rivi
    line2 <- fid[2] # toinen rivi
    line3 <- fid[3] # kolmas
  } else {
    fid <- line1 <- line2 <- line3 <- -1
  }

  if (line1 == -1 | line2 == -1 | line3 == -1) {
    stop("Incorrect file format 1168")
  }
  if (testaaPop(line1) == 1 | testaaPop(line2) == 1) {
    stop("Incorrect file format 1172")
  }
  if (testaaPop(line3) == 1) {
    # 2 rivi t�ll�in lokusrivi (2 rows then locus row)
    nloci <- rivinSisaltamienMjonojenLkm(line2)
    line4 <- fid[4]
    if (line4 == -1) stop("Incorrect file format 1180")
    if (!grepl(",", line4)) {
      # Rivin nelj?t�ytyy sis�lt�� pilkku.
      stop("Incorrect file format 1185")
    }
    pointer <- 1
    while (substring(line4, pointer, pointer) != ",") {
      # Tiedet��n, ett?pys�htyy
      pointer <- pointer + 1
    }

    # pilkun j�lkeinen osa (the part after the comma)
    line4 <- substring(line4, pointer + 1)

    nloci2 <- rivinSisaltamienMjonojenLkm(line4)
    if (nloci2 != nloci) stop("Incorrect file format 1195")
  } else {
    line <- fid[4]
    lineNumb <- 4
    while (testaaPop(line) != 1 & line != -1) {
      line <- fid[lineNumb + 1]
      lineNumb <- lineNumb + 1
    }
    if (line == -1) stop("Incorrect file format 1206")
    nloci <- lineNumb - 2
    line4 <- fid[lineNumb + 1] # Eka rivi pop sanan j�lkeen
    if (line4 == -1) stop("Incorrect file format 1212")
    if (!grepl(",", line4)) {
      # Rivin t�ytyy sis�lt�� pilkku. (The line must contain a comma)
      stop("Incorrect file format 1217")
    }
    pointer <- 1
    while (substring(line4, pointer, pointer) != ",") {
      # Tiedet��n, ett?pys�htyy
      pointer <- pointer + 1
    }

    # pilkun j�lkeinen osa (the part after the comma)
    line4 <- substring(line4, pointer + 1)

    nloci2 <- rivinSisaltamienMjonojenLkm(line4)
    if (nloci2 != nloci) stop("Incorrect file format 1228")
  }
  kunnossa <- 1
  return(kunnossa)
}
