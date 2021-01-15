#' @title suoritaMuutos
#' @description Päivittää osuusTaulun muutoksen jälkeen.
#' @param osuusTaulu Percentage table?
#' @param osuus percentage?
#' @param indeksi index
#' @export
suoritaMuutos <- function (osuusTaulu, osuus, indeksi) {
	if (is.null(dim(COUNTS))) {
		npops <- 1
	} else {
		npops <- ifelse(is.na(dim(COUNTS)[3]), 1, dim(COUNTS)[3])
	}

	i1 <- indeksi %% npops
	if (is.na(i1) | i1 == 0) i1 <- npops
	i2 <- ceiling(indeksi / npops)

	osuusTaulu[i1] <- osuusTaulu[i1] - osuus
	osuusTaulu[i2] <- osuusTaulu[i2] + osuus

	return (osuusTaulu)
}