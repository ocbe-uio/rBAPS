#' @title Take line
#' @description Returns one line from the description.
#' @param description description
#' @param width width
#' @return newline
takeLine <- function(description, width) {
  # Returns one line from the description: line ends to the first
  # space after width:th mark.
  n <- width + 1
  while (description[n] != "" && n < length(description)) {
    n <- n + 1
  }
  newline <- description[1:n]
  return(newline)
}
