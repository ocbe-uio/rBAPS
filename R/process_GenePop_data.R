process_GenePop_data <- function(filename) {
  # Pre-processing GenePop data
  kunnossa <- testaaGenePopData(filename)
  data_popnames <- lueGenePopData(filename)
  data <- data_popnames$data
  popnames <- data_popnames$popnames

  # Processing GenePop data
  c <- handleData(data, "genepop")
  Z_dist <- newGetDistances(c$newData, c$rowsFromInd)

  # Forming and returning pre-processed data
  list(
    "data" = c$newData,
    "rowsFromInd" = c$rowsFromInd,
    "alleleCodes" = c$alleleCodes,
    "noalle" = c$noalle,
    "adjprior" = c$adjprior,
    "priorTerm" = c$priorTerm,
    "dist" = Z_dist$dist,
    "popnames" = popnames,
    "Z" = Z_dist$Z
  )
}
