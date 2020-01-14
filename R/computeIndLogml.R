#' @title computeIndLogml
#' @description Palauttaa yksilön logml:n, kun oletetaan yksilön alkuperät
#' määritellyiksi kuten osuusTaulu:ssa.
#' @param omaFreqs omaFreqs
#' @param osuusTaulu osuusTaulu
#' @export
computeIndLogml <- function (omaFreqs, osuusTaulu) {

    apu <- repmat(t(osuusTaulu), c(1, dim(omaFreqs)[2]))
    apu <- c(apu) * omaFreqs # c() avoids deprecation error re. matrix ops
    if (length(apu) > 1) {
        apu <- colSums(as.matrix(apu))
    } else {
        apu <- sum(apu)
    }

    apu = log(apu)

    loggis <- sum(apu)
    return (loggis)
}