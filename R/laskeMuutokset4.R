#' @title laskeMuutokset4
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
    notEmpty <- osuusTaulu > 0.005
    muutokset <- zeros(npops)
    empties <- !notEmpty

    for (i1 in notEmpty) {
        osuusTaulu[i1] <- osuusTaulu[i1] - osuus
        for (i2 in c(1:(i1 - 1), (i1 + 1):npops)) {
            osuusTaulu[i2] <- osuusTaulu[i2] + osuus
            loggis <- computeIndLogml(omaFreqs, osuusTaulu)
            muutokset[i1, i2] <- loggis - logml
            osuusTaulu[i2] <- osuusTaulu[i2] - osuus
        } 
        osuusTaulu[i1] <- osuusTaulu[i1] + osuus
    }
    return (muutokset)
}