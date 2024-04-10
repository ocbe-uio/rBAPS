computePopulationLogml <- function(pops, adjprior, priorTerm) {
  # Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset

  nr <- seq_len(nrow(adjprior))
  nc <- seq_len(ncol(adjprior))

  x <- size(baps.globals$COUNTS, 1)
  y <- size(baps.globals$COUNTS, 2)
  z <- length(pops)

  # ======================================================== #
  # Computation                                              #
  # ======================================================== #
  rep_adj <- repmat(adjprior, c(1, 1, z))
  gamma_rep_counts <- matlab2r::gammaln(rep_adj + baps.globals$COUNTS[, , pops])
  gamma_sum_counts <- rowSums(matlab2r::gammaln(1 + baps.globals$SUMCOUNTS[pops, , drop = FALSE]))
  gamma_rep_counts_sum <- colSums(colSums(reshape(gamma_rep_counts, c(x, y, z))))
  gamma_rep_counts_reshaped <- squeeze(gamma_rep_counts_sum)
  popLogml <- gamma_rep_counts_reshaped - gamma_sum_counts - priorTerm
  return(popLogml[, , drop = FALSE])
}
