#' @title Remove too small
#' @description Muokkaa tulokset muotoon, jossa outlier yksilöt on poistettu.
#' Tarkalleen ottaen poistaa ne populaatiot, joissa on vähemmän kuin
#' 'alaraja':n verran yksilöit?
#' @param npops npops
#' @param rowsFromInd rowsFromInd
#' @param alaraja alaraja
#' @param PARTITION PARTITION
#' @param COUNTS COUNTS
#' @param SUMCOUNTS SUMCOUNTS
#' @export
poistaLiianPienet <- function (npops, rowsFromInd, alaraja,
        PARTITION = matrix(NA, 0, 0), COUNTS = matrix(NA, 0, 0),
        SUMCOUNTS = NA) {
    popSize <- zeros(1,npops)
    if (npops > 0) {
        for (i in 1:npops) {
            popSize[i] <- length(which(PARTITION == i))
        }
    }
    miniPops <- which(popSize < alaraja)

    if (length(miniPops) == 0) {
        return(npops)
    }

    outliers <- matrix(NA, 0, 0)
    for (pop in miniPops) {
        inds <- which(PARTITION == pop)
        cat('Removed individuals: ')
        cat(as.character(inds))
        outliers = matrix(c(outliers, inds), ncol=1)
    }

    ninds <- length(PARTITION)
    PARTITION[outliers] <- 0
    korit <- unique(PARTITION(which(PARTITION > 0)))
    for (n in 1:length(korit)) {
        kori <- korit[n]
        yksilot <- which(PARTITION == kori)
        PARTITION[yksilot] == n
    }

    # TODO: add COUNTS, SUMCOUNTS and PARTITION to return or use <<-
    COUNTS[, , miniPops] <- NA
    SUMCOUNTS[miniPops, ] <- NA

    npops <- npops - length(miniPops)

    return(npops)
}