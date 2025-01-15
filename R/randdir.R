#' @title Generates random numbers
#' @return vector of length `nc` with r.v. realizations from Gamma(rate=1)
#' @examples rBAPS:::randdir(matrix(c(10, 30, 60), 3), 3)
#' @param counts shape parameter
#' @param nc number of rows on output
#' @importFrom stats rgamma
randdir <- function(counts, nc) {
  svar <- zeros(nc, 1)
  for (i in 1:nc) {
    svar[i, 1] <- rgamma(1L, shape = counts[i, 1], rate = 1)
  }
  svar <- svar / sum(svar)
  return(svar)
}
