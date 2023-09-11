poistaTyhjatPopulaatiot <- function(npops) {
  # % Poistaa tyhjentyneet populaatiot COUNTS:ista ja
  # % SUMCOUNTS:ista. P�ivitt�� npops:in ja PARTITION:in.
  notEmpty <- matlab2r::find(apply(SUMCOUNTS, 1, function(x) any(x > 0)))
  COUNTS <- COUNTS[, , notEmpty]
  SUMCOUNTS <- SUMCOUNTS[notEmpty, ]
  LOGDIFF <- LOGDIFF[, notEmpty]

  for (n in 1:length(notEmpty)) {
    apu <- matlab2r::find(PARTITION == notEmpty[n])
    PARTITION[apu] <- n
  }
  npops <- length(notEmpty)
  return(npops)
}
