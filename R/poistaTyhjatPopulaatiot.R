poistaTyhjatPopulaatiot <- function(npops) {
  # % Poistaa tyhjentyneet populaatiot COUNTS:ista ja
  # % SUMCOUNTS:ista. P�ivitt�� npops:in ja PARTITION:in.
  notEmpty <- find(any(SUMCOUNTS, 2))
  COUNTS <- COUNTS[, , notEmpty]
  SUMCOUNTS <- SUMCOUNTS[notEmpty, ]
  LOGDIFF <- LOGDIFF[, notEmpty]

  for (n in 1:length(notEmpty)) {
    apu <- find(PARTITION == notEmpty(n))
    PARTITION[apu] <- n
  }
  npops <- length(notEmpty)
  return(npops)
}
