#' @title Logml to string
#' @description Returns a string representation of a logml
#' @param logml input Logml
#' @param
#' @param leading_zeros_replacement string to replace leading zeros with
#' @return String version of logml
logml2String <- function(logml, leading_zeros_replacement = " ") {
  mjono <- rep(" ", 7L)
  # Palauttaa logml:n string-esityksen.

  if (logml == -Inf) {
    mjono[7] <- "-"
    return(mjono)
  }

  if (abs(logml) < 10000) {
    # Ei tarvita e-muotoa
    mjono[7] <- palautaYks(abs(logml), -1)
    mjono[6] <- "."
    mjono[5] <- palautaYks(abs(logml), 0)
    mjono[4] <- palautaYks(abs(logml), 1)
    mjono[3] <- palautaYks(abs(logml), 2)
    mjono[2] <- palautaYks(abs(logml), 3)
    pointer <- 2
    while (mjono[pointer] == "0" && pointer < 7) {
      # Removes leading zeros
      mjono[pointer] <- leading_zeros_replacement
      pointer <- pointer + 1
    }
    if (logml < 0) {
      mjono[pointer - 1] <- "-"
    }
  } else {
    suurinYks <- 4
    while (abs(logml) / (10^(suurinYks + 1)) >= 1) {
      suurinYks <- suurinYks + 1
    }
    if (suurinYks < 10) {
      mjono[7] <- as.character(suurinYks)
      mjono[6] <- "e"
      mjono[5] <- palautaYks(abs(logml), suurinYks - 1)
      mjono[4] <- "."
      mjono[3] <- palautaYks(abs(logml), suurinYks)
      if (logml < 0) {
        mjono[2] <- "-"
      }
    } else if (suurinYks >= 10) {
      mjono[6:7] <- as.character(suurinYks)
      mjono[5] <- "e"
      mjono[4] <- palautaYks(abs(logml), suurinYks - 1)
      mjono[3] <- "."
      mjono[2] <- palautaYks(abs(logml), suurinYks)
      if (logml < 0) {
        mjono[1] <- "-"
      }
    }
  }
  return(mjono)
}
