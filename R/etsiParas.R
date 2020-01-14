etsiParas <- function = (osuus, osuusTaulu, omaFreqs, logml) {
    ready <- 0;
    while (ready != 1) {
        muutokset <- laskeMuutokset4(osuus, osuusTaulu, omaFreqs, logml)
        [maxMuutos, indeksi] = max(muutokset[1:end]) # TODO: how does this work on Matlab?
        if (maxMuutos > 0) {
            osuusTaulu <- suoritaMuutos(osuusTaulu, osuus, indeksi)
            logml <- logml + maxMuutos
        } else {
            ready <- 1
        }
    }
    return (c(osuusTaulu, logml))
}