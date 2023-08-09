#' @title simuloiAlleeli
#' @description Simuloi populaation pop lokukseen loc alleelin.
#' @note This function is (only?) called by `simulateIndividuals()`. Therefore,
#' exporting it is probably unnecessary.
#' @param allfreqs allfreqa
#' @param pop pop
#' @param loc loc
simuloiAlleeli <- function(allfreqs, pop, loc) {
  if (length(dim(allfreqs)) == 0) {
    freqs <- 1
  } else {
    if (length(dim(allfreqs)) == 3) { # distinguish between array and matrix
      freqs <- allfreqs[, loc, pop]
    } else {
      freqs <- allfreqs[, loc]
    }
  }
  # freqs <- ifelse(is.null(length(dim(allfreqs)), allfreqs[loc], 0)
  # freqs <- switch() + 1,
  #     allfreqs[, loc],
  #     allfreqs[, loc, pop]
  # )


  cumsumma <- cumsum(freqs)
  arvo <- runif(1)
  isommat <- which(cumsumma > arvo)
  all <- base::min(isommat)
  return(all)
}
