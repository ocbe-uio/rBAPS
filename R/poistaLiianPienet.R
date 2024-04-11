#' @title Remove too small
#' @description Muokkaa tulokset muotoon, jossa outlier yksilöt on poistettu.
#' Tarkalleen ottaen poistaa ne populaatiot, joissa on vähemmän kuin
#' 'alaraja':n verran yksilöit?
#' @param npops npops
#' @param rowsFromInd rowsFromInd
#' @param alaraja alaraja
poistaLiianPienet <- function(npops, rowsFromInd, alaraja) {
  popSize <- zeros(1, npops)
  if (npops > 0) {
    for (i in 1:npops) {
      popSize[i] <- length(which(globals$PARTITION == i))
    }
  }
  miniPops <- which(popSize < alaraja)

  if (length(miniPops) == 0) {
    return(npops)
  }

  outliers <- matrix(NA, 0, 0)
  for (pop in miniPops) {
    inds <- which(globals$PARTITION == pop)
    cat("Removed individuals: ")
    cat(as.character(inds))
    outliers <- matrix(c(outliers, inds), ncol = 1)
  }

  ninds <- length(globals$PARTITION)
  globals$PARTITION[outliers] <- 0
  korit <- unique(globals$PARTITION(which(globals$PARTITION > 0)))
  for (n in 1:length(korit)) {
    kori <- korit[n]
    yksilot <- which(globals$PARTITION == kori)
    globals$PARTITION[yksilot] == n
  }

  # TODO: add COUNTS, SUMCOUNTS and PARTITION to return or use <-
  COUNTS[, , miniPops] <- NA
  SUMCOUNTS[miniPops, ] <- NA

  npops <- npops - length(miniPops)

  return(npops)
}
