fiksaaPartitioYksiloTasolle <- function(rows, rowsFromInd) {
  # Fix partition to individual level
  totalRows <- 0
  for (ind in 1:size(rows, 1)) {
    totalRows <- totalRows + (rows[ind, 2] - rows[ind, 1] + 1)
  }
  partitio2 <- zeros(totalRows / rowsFromInd, 1)

  for (ind in 1:size(rows, 1)) {
    kaikkiRivit <- rows[ind, 1]:rows[ind, 2]
    for (riviNumero in seq(rowsFromInd, length(kaikkiRivit), rowsFromInd)) {
      rivi <- kaikkiRivit[riviNumero]
      partitio2[rivi / rowsFromInd] <- PARTITION[ind]
    }
  }
  global_env <- as.environment(1L)
  assign("PARTITION", partitio2, envir = global_env)
}
