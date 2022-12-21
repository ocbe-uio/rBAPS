laskeKlikit <- function(M, maxCliqSize, maxSepSize) {
  # Laskee samankokoisten klikkien mההrהn verkosta M
  # ncliques(i)=kokoa i olevien klikkien mההr?
  # nseparators vastaavasti
  ncliques <- zeros(1, maxCliqSize)
  nseparators <- zeros(1, maxSepSize)
  if (M == c()) {
    return()
  }
  cliques_separators <- findCliques(M)
  cliques <- cliques_separators$cliques
  separators <- cliques_separators$separators
  rm(cliques_separators)
  for (i in 1:length(cliques)) {
    ncliques[length[cliques[[i]]]] <- ncliques[length(cliques[[i]])] + 1
  }
  # cliqmax=max(find(ncliques!=0))
  # ncliques=ncliques(1:cliqmax)
  for (i in 1:length(separators)) {
    nseparators[length[separators[[i]]]] <- nseparators[length(separators[[i]])] + 1
  }
  return(
    list(
      ncliques = ncliques, nseparators = nseparators, cliques = cliques,
      separators = separators
    )
  )
}
