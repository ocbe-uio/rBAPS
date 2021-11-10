initializeGammaln <- function(ninds, rowsFromInd, maxAlleles) {
  # Alustaa GAMMALN muuttujan s.e. GAMMALN(i, j)=gammaln((i - 1) + 1/j)
  GAMMA_LN <- zeros((1 + ninds) * rowsFromInd, maxAlleles)
  for (i in 1:(ninds + 1) * rowsFromInd) {
    for (j in 1:maxAlleles) {
      GAMMA_LN[i, j] <- log_gamma((i - 1) + 1 / j)
    }
  }
}
