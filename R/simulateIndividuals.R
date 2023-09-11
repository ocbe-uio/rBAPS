#' @title Simulate individuals
#' @description simulate n individuals from population pop, such that
#' proportion "missing_level" of the alleles are present.
#' @param n n
#' @param rowsFromInd rowsFromInd
#' @param allfreqs allfreqs
#' @param pop pop
#' @param missing_level missing_level
simulateIndividuals <- function(n, rowsFromInd, allfreqs, pop, missing_level) {
  nloci <- size(allfreqs, 2)

  refData <- zeros(n * rowsFromInd, nloci)
  counter <- 1 # which row will be generated next.

  for (ind in 1:n) {
    for (loc in 1:nloci) {
      for (k in 0:(rowsFromInd - 1)) {
        if (runif(1) < missing_level) {
          refData[counter + k, loc] <- simuloiAlleeli(
            allfreqs, pop, loc
          )
        } else {
          refData[counter + k, loc] <- -999
        }
      }
    }
    counter <- counter + rowsFromInd
  }
  return(refData)
}
