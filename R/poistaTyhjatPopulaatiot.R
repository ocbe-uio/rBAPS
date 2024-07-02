poistaTyhjatPopulaatiot <- function(npops) {
  # % Poistaa tyhjentyneet populaatiot globals$COUNTS:ista ja
  # % globals$SUMCOUNTS:ista. P�ivitt�� npops:in ja globals$PARTITION:in.
  notEmpty <- matlab2r::find(apply(globals$SUMCOUNTS, 1, function(x) any(x > 0)))
  globals$COUNTS <- globals$COUNTS[, , notEmpty]
  globals$SUMCOUNTS <- globals$SUMCOUNTS[notEmpty, ]
  globals$LOGDIFF <- globals$LOGDIFF[, notEmpty]

  for (n in 1:length(notEmpty)) {
    apu <- matlab2r::find(globals$PARTITION == notEmpty[n])
    globals$PARTITION[apu] <- n
  }
  npops <- length(notEmpty)
  return(npops)
}
