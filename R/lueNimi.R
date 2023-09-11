#' @title Read the Name
#' @description Returns the part of the line from the beginning that is before the comma. Useful for returning the name of a GenePop area
#' @param line line
#' @return nimi
lueNimi <- function(line) {
  # ==========================================================================
  # Validation
  # ==========================================================================
  if (!grepl(",", line)) {
    stop("There are no commas in this line")
  }
  # Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
  n <- 1
  merkki <- substring(line, n, n)
  nimi <- ""
  while (merkki != ",") {
    nimi <- c(nimi, merkki)
    n <- n + 1
    merkki <- substring(line, n, n)
  }
  return(paste(nimi, collapse = ""))
}
