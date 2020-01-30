#' @title Compute Personal Freqs
#' @description Laskee npops*(rowsFromInd*nloci) taulukon, jonka kutakin
#' saraketta vastaa yksilön ind alleeli. Eri rivit ovat alleelin
#' alkuperäfrekvenssit eri populaatioissa. Jos yksilölt?puuttuu jokin alleeli, 
#' niin vastaavaan kohtaa tulee sarake ykkösi?
#' @param ind ind
#' @param data data
#' @param allFreqs allFreqs
#' @param rowsFromInd rowsFromInd
#' @param COUNTS COUNTS
#' @export

computePersonalAllFreqs <- function(ind, data, allFreqs, rowsFromInd,
COUNTS = matrix(0)) {
    nloci <- ifelse(is.na(dim(COUNTS)[2]), 1, dim(COUNTS)[2])
    npops <- ifelse(is.na(dim(COUNTS)[3]), 1, dim(COUNTS)[3])

    rows <- as.matrix(t(data))[computeRows(rowsFromInd, ind, 1), , drop = FALSE]

    omaFreqs <- zeros(npops, rowsFromInd * nloci)
    pointer <- 1
    for (loc in 1:dim(rows)[2]) {
        for (all in 1:dim(rows)[1]) {
            if (rows[all, loc] >= 0) {
                if (pointer > ncol(omaFreqs)) omaFreqs <- cbind(omaFreqs, 0)
                omaFreqs[, pointer] <- tryCatch(
                    matrix(
                        data = as.matrix(t(allFreqs))[rows[all, loc], loc],
                        nrow = npops
                    ),
                    error = function(e) return(NA)
                )
            } else {
                omaFreqs[, pointer] <- ones(npops, 1)
            }
            # omaFreqs <- unname(cbind(omaFreqs, new_omaFreqs))
            pointer <- pointer + 1
        }
    }
    omaFreqs <- omaFreqs[, !is.na(omaFreqs)]
    return(omaFreqs)
}