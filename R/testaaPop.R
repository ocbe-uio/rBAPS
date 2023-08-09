#' @title Test population
#' @description Test a line in the population
#' @param rivi Line
#' @return pal = 1 if the line starts with one of the following
# letter combinations: Pop, pop, POP. In all others cases, pal = 0
testaaPop <- function(rivi) {
  # pal=1, mik�li rivi alkaa jollain seuraavista
  # kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
  # tapauksissa pal=0.

  if (nchar(rivi) < 3) {
    pal <- 0
  } else {
    rivi_start <- substring(rivi, 1, 3)
    pal <- ifelse(rivi_start %in% c("Pop", "pop", "POP"), 1, 0)
  }
  return(pal)
}
