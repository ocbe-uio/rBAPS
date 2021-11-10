#' @title computeIndLogml
#' @description Palauttaa yksilön logml:n, kun oletetaan yksilön alkuperät
#' määritellyiksi kuten osuusTaulu:ssa.
#' @param omaFreqs own Freqs?
#' @param osuusTaulu Percentage table?
#' @export
computeIndLogml <- function(omaFreqs, osuusTaulu) {
  omaFreqs <- as.matrix(omaFreqs)
  osuusTaulu <- as.matrix(osuusTaulu)

  apu <- repmat(t(osuusTaulu), c(1, dim(omaFreqs)[2]))
  apu <- times(apu, omaFreqs) # c() avoids deprecation error re. matrix ops
  if (length(apu) > 1) {
    apu <- colSums(as.matrix(apu))
  } else {
    apu <- sum(apu)
  }

  if (any(apu < 0)) {
    # Workaround for log of a negative number
    apu <- as.complex(apu)
  }
  apu <- log(apu)

  loggis <- sum(apu)
  return(loggis)
}
