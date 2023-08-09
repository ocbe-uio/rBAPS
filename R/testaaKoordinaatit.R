#' @title Test the coordinates
#' @param ninds ninds
#' @param coordinates coordinates
#' @param interactive prompt user for relevant questions during execution
#' @return a list of defectives ("viallinen") and coordinates
testaaKoordinaatit <- function(ninds, coordinates, interactive = TRUE) {
 # Testaa onko koordinaatit kunnollisia.
 # modified by Lu Cheng, 05.12.2012
  viallinen <- 1
  if (any(sapply(coordinates, class) != "numeric")) {
    warning('Coordinates are not numerical!')
    return()
  }
  oikeanKokoinen <- size(coordinates, 1) == ninds & (size(coordinates, 2) == 2)
  if (!oikeanKokoinen) {
    warning('Wrong coordinates dimension!')
    return()
  }
  posstr <- sapply(coordinates, function(x) sprintf('%.10f', x))
  posstr <- gsub('\\.0.+', '.', posstr)
  posstr <- matrix(posstr, nrow = nrow(coordinates))
  uni1 <- unique(posstr[, 1])
  uni2 <- unique(posstr[, 2])
  posstr_new <- posstr
  if (length(uni1) == ninds && length(uni2) == ninds) {
    viallinen <- 0
    return(list(viallinen = viallinen, coordinates = coordinates))
  } else {
    ans <- "Yes"
    if (interactive) {
      ans <- questdlg(
        'Input coordinates are not unique. Do you want to make them unique?',
        'coordinates NOT unique', c('Yes', 'No'), 'Yes'
      )
    }
    if (strcmp(tolower(ans), 'no')) {
      warning('Coordinates are not unique!')
      return(list(viallinen = viallinen, coordinates = coordinates))
    }
  }

  for (i in 1:length(uni1)) {
    tmpinds <- find(posstr[, 1] %in% uni1[i])
    tmpNinds <- length(tmpinds)
    if (tmpNinds == 1) {
      next
    }
    if (tmpNinds >= 100) stop("Assertion failed. tmpNinds not < 100")
    tmparr <- round(seq(0, 99, length.out = tmpNinds))
    tmparr <- tmparr[sample(tmpNinds)]
    for (j in 1:tmpNinds) {
      posstr_new[tmpinds[j], 1] <- sprintf('%s%02d', posstr[tmpinds[j], 1], tmparr[j])
    }
  }

  for (i in 1:length(uni2)) {
    tmpinds <- find(posstr[, 2] %in% uni2[i])
    tmpNinds <- length(tmpinds)
    if (tmpNinds == 1) next
    if (tmpNinds >= 100) stop("Assertion failed. tmpNinds not < 100")
    tmparr <- round(seq(0, 99, length.out = tmpNinds))
    tmparr <- tmparr[sample(tmpNinds)]
      for (j in 1:tmpNinds) {
          posstr_new[tmpinds[j], 2] <- sprintf('%s%02d', posstr[tmpinds[j], 2], tmparr[j])
      }
  }
  coordinates <- matrix(sapply(posstr_new, as.double), ncol = 2)
  uni1 <- unique(coordinates[, 1])
  uni2 <- unique(coordinates[, 2])
  if (length(uni1 )== ninds && length(uni2) == ninds) {
      viallinen <- 0
  } else {
      warning('Can not make coordinates unique!')
  }
  return(list(viallinen = viallinen, coordinates = coordinates))
}
