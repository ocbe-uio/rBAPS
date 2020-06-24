#' @title Test population
#' @description Test a line in the population
#' @param rivi Line
#' @return pal
#' @export
testaaPop <- function(rivi) {
	# pal=1, mik�li rivi alkaa jollain seuraavista
	# kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
	# tapauksissa pal=0.

	if (length(rivi) < 3) {
		pal <- 0
		return(pal)
	}
	if (
		all(rivi[1:3] == 'Pop') |
		all(rivi[1:3] == 'pop') |
		all(rivi(1:3)=='POP')
	) {
		pal <- 1
		return(pal)
	} else {
		pal <- 0;
		return(pal)
	}
	return(pal)
}