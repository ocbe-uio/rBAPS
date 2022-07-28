computeDiffInCliqCounts <- function(cliques, inds) {
  # Laskee muutoksen CLIQCOUNTS:ssa (tai SEPCOUNTS:ssa, jos syצtteen?
  # separators) kun yksilצt inds siirretההn. diffInCliqcounts on ncliq*1 taulu,
  # joka on CLIQCOUNTS:n sarakkeesta josta yksilצt inds siirretההn ja
  # lisהttהv?sarakkeeseen, johon yksilצt siirretההn.
  ncliq <- size(cliques, 1)
  diffInCliqCounts <- zeros(ncliq, 1)
  ninds <- length(inds)
  for (i in 1:ninds) {
    ind <- inds[i]
    rivit <- rowSums(cliques == ind)
    diffInCliqCounts <- diffInCliqCounts + rivit
  }
  return(diffInCliqCounts)
}
