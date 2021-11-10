kldiv2str <- function(div) {
  mjono <- "      "
  if (abs(div) < 100) {
    # Ei tarvita e-muotoa
    mjono[6] <- as.character((floor(div * 1000)) %% 10)
    mjono[5] <- as.character((floor(div * 100)) %% 10)
    mjono[4] <- as.character((floor(div * 10)) %% 10)
    mjono[3] <- "."
    mjono[2] <- as.character((floor(div)) %% 10)
    arvo <- (floor(div / 10)) %% 10
    if (arvo > 0) {
      mjono[1] <- as.character(arvo)
    }
  } else {
    suurinYks <- floor(log10(div))
    mjono[6] <- as.character(suurinYks)
    mjono[5] <- "e"
    mjono[4] <- palautaYks(abs(div), suurinYks - 1)
    mjono[3] <- "."
    mjono[2] <- palautaYks(abs(div), suurinYks)
  }
  return(mjono)
}
