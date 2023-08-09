#' @title Initialize Pop Names
#' @param nameFile nameFile
#' @param indexFile indexFile
initPopNames <- function(nameFile, indexFile) {
  # Palauttaa tyhj�n, mik�li nimitiedosto ja indeksitiedosto
  #  eiv�t olleet yht?pitki?

  indices <- load(indexFile)

  fid <- load(nameFile)
  if (fid == -1) {
    # File didn't exist
    stop("Loading of the population names was unsuccessful")
  }
  line <- readLines(fid)[1]
  counter <- 1
  names <- vector()
  while ((line != -1) & (line != "")) {
    names[counter] <- line
    line <- readLines(fid)[counter]
    counter <- counter + 1
  }

  if (length(names) != length(indices)) {
    cat("The number of population names must be equal to the number of")
    cat("entries in the file specifying indices of the first individuals")
    cat("of each population.")
  }

  popnames <- cell(length(names), 2)
  for (i in 1:length(names)) {
    popnames[i, 1] <- names[i]
    popnames[i, 2] <- indices[i]
  }
  return(popnames)
}
