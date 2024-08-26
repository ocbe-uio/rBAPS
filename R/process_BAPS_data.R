process_BAPS_data <- function(file, partitionCompare) {
  if (!is.null(partitionCompare)) {
    cat('Data:', file, '\n')
  }

  # Importing data
  if (is.character(file)) {
    data <- read.table(file)
  } else {
    data <- file
  }

  ninds <- testaaOnkoKunnollinenBapsData(data)  # Checks if last column is ID
  if (ninds == 0) {
    warning('Incorrect Data-file.')
    return(NULL)
  }

  popnames <- NULL  # Dropped specification of population names (from BAPS 6)

  # Processing data
  result <- handleData(data, format = "BAPS")
  data <- result$newData
  rowsFromInd <- result$rowsFromInd
  alleleCodes <- result$alleleCodes
  noalle <- result$noalle
  adjprior <- result$adjprior
  priorTerm <- result$priorTerm
  result <- newGetDistances(data, rowsFromInd)
  Z <- result$Z
  dist <- result$dist

  # Forming and saving pre-processed data
  processed_data <- list(
    data = data,
    rowsFromInd = rowsFromInd,
    alleleCodes = alleleCodes,
    noalle = noalle,
    adjprior = adjprior,
    priorTerm = priorTerm,
    dist = dist,
    popnames = popnames,
    Z = Z
  )
  return(processed_data)
}
