#' @title simuloiAlleeli
#' @description Simuloi populaation pop lokukseen loc alleelin.
#' @note This function is (only?) called by `simulateIndividuals()`. Therefore, exporting it is probably unnecessary.
#' @export

simuloiAlleeli <- function(allfreqs, pop, loc) {
    if (length(dim(allfreqs)) == 3) { # distinguish between arrays and matrices
        freqs <- allfreqs[, loc, pop]
    } else {
        freqs <- allfreqs[, loc]
    }
    cumsumma <- cumsum(freqs)
    arvo <- runif(1)
    isommat <- which(cumsumma > arvo)
    all <- min(isommat)
    return(all)
}