#' @title Calculate changes?
#' @description Palauttaa npops*npops taulun, jonka alkio (i,j) kertoo, mik?on
#' muutos logml:ss? mikäli populaatiosta i siirretään osuuden verran
#' todennäköisyysmassaa populaatioon j. Mikäli populaatiossa i ei ole mitään 
#' siirrettävää, on vastaavassa kohdassa rivi nollia.
#' @param osuus Percentages?
#' @param omaFreqs own Freqs?
#' @param osuusTaulu Percentage table?
#' @param logml log maximum likelihood
#' @param COUNTS COUNTS
#' @export
laskeMuutokset4 <- function (osuus, osuusTaulu, omaFreqs, logml,
                             COUNTS = matrix(0)) {
    npops <- ifelse(is.na(dim(COUNTS)[3]), 1, dim(COUNTS)[3])
    notEmpty <- which(osuusTaulu > 0.005)
    muutokset <- zeros(npops)
    empties <- !notEmpty

    for (i1 in notEmpty) {
        osuusTaulu[i1] <- osuusTaulu[i1] - osuus
        for (i2 in c(colon(1, i1 - 1), colon(i1 + 1, npops))) {
            osuusTaulu[i2] <- osuusTaulu[i2] + osuus
            loggis <- computeIndLogml(omaFreqs, osuusTaulu)
            
            # Work around Matlab OOB bug
            if (i1 > nrow(muutokset)) {
                muutokset <- rbind(muutokset, muutokset * 0)
            }
            if (i2 > ncol(muutokset)) {
                muutokset <- cbind(muutokset, muutokset * 0)
            }

            muutokset[i1, i2] <- loggis - logml
            osuusTaulu[i2] <- osuusTaulu[i2] - osuus
        } 
        osuusTaulu[i1] <- osuusTaulu[i1] + osuus
    }
    return (muutokset)
}