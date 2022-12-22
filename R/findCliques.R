findCliques <- function(M) {
  # Muuttaa graafin M kolmioituvaksi ja laskee siitה klikit ja
  # separaattorit.
  # Hyצdynnetההn Kevin Murphyn algoritmeja Graph Theory toolboxista.
  # Pהivitetty 12.8.2005
  order <- elim_order(M, ones(length(M))) # TODO: translate from findCliques.m
  G <- cliques <- root <- NULL
  c(G, cliques) %<-% triangulate(M, order) # TODO: translate from findCliques.m
  c(jtree, root) %<-% cliques_to_jtree(cliques, ones(length(M))) # TODO: translate from findCliques.m
  ncliq <- length(cliques)
  separators <- cell(ncliq - 1, 1)    # n - solmuisessa puussa n - 1 viivaa

  jono <- zeros(length(ncliq))
  jono[1] <- root
  i <- 1
  pointer <- 2                     # Seuraava tyhjה paikka

  while (!is.null(find(jono != 0))) {  # Puun leveyssuuntainen lהpikהynti) {
    lapset <- find(jtree[jono[i], ] != 0)
    jtree[, jono[i]] <- 0         # Klikki kהsitelty
    jono[pointer:(pointer + length(lapset) - 1)] <- lapset
    for (j in 1:length(lapset)) {
      ehdokas <- myintersect(cliques[[jono[i]]], cliques[[lapset[j]]])
      kelpaa <- 1
      for (k in 1:(pointer + j - 3)) {
        # Tutkitaan, ettה separaattoriehdokasta ei vielה kהsitelty
        if (ehdokas == separators[[k]]) {
          kelpaa <- 0
        }
      }
      if (kelpaa) {
        separators[[pointer + j - 2]] <- ehdokas
      }
    }
    jono[i] <- 0
    pointer <- pointer + length(lapset)
    i <- i + 1
  }
  notEmpty <- zeros(ncliq - 1, 1)
  for (i in 1:(ncliq - 1)) {
    if (!is.null(separators[[i]])) {
      notEmpty[i] <- 1
    }
  }
  notEmpty <- find(notEmpty == 1)
  separators <- separators(notEmpty)
  return(list("cliques" = cliques, "separators" = separators, "G" = G))
}
