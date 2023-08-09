computePopulationLogml <- function(pops, adjprior, priorTerm) {
  # Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset

  # ======================================================== #
  # Limiting COUNTS size                                     #
  # ======================================================== #
  if (!is.null(adjprior)) {
    nr <- seq_len(nrow(adjprior))
    nc <- seq_len(ncol(adjprior))
    COUNTS <- COUNTS[nr, nc, pops, drop = FALSE]
  } else {
    COUNTS <- NA
  }

  x <- size(COUNTS, 1)
  y <- size(COUNTS, 2)
  z <- length(pops)

  # ======================================================== #
  # Computation                                              #
  # ======================================================== #
  term1 <- NULL
  if (!is.null(adjprior)) {
    isarray <- length(dim(repmat(adjprior, c(1, 1, length(pops))))) > 2
    term1 <- squeeze(
      sum(
        sum(
          reshape(
            lgamma(
              repmat(adjprior, c(1, 1, length(pops))) + COUNTS[nr, nc, pops, drop = !isarray]
            ),
            c(x, y, z)
          ),
          1
        ),
        2
      )
    )
  }
  if (is.null(priorTerm)) priorTerm <- 0
  popLogml <- term1 - sum(lgamma(1 + SUMCOUNTS[pops, ]), 2) - priorTerm
  return(popLogml)
}
