laskeVarmuus <- function(
  rowsFromInd, data, adjprior, priorTerm, logml, cliques, separators, ninds
) {
  varmuus <- zeros(ninds, 1)
  for (ind in 1:ninds) {
    muutokset <- spatialMixture_muutokset$new()
    muutokset <- muutokset$laskeMuutokset(
      ind, rowsFromInd, data, adjprior, priorTerm, logml, cliques, separators
    )
    varmuus[ind] <- 1 / sum(exp(muutokset))
  }
  return(varmuus)
}
