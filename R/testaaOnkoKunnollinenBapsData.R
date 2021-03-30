#' @title TestBAPS data
#' @description Test if loaded BAPS data is proper
#' @param data dataset
#' @return ninds
#' @export
testaaOnkoKunnollinenBapsData <- function(data) {
	# Tarkastaa onko viimeisess?sarakkeessa kaikki
	# luvut 1,2,...,n johonkin n:��n asti.
	# Tarkastaa lis�ksi, ett?on v�hint��n 2 saraketta.
	if (size(data, 1) < 2) {
		ninds <- 0
		return(ninds)
	}
	lastCol <- data[, ncol(data)]
	ninds <- max(lastCol)
	if (any(1:ninds != unique(lastCol))) {
		ninds <- 0
		return(ninds)
	}
	return(ninds)
}