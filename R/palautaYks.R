palautaYks <- function(num, yks) {
	# palauttaa luvun num 10^yks termin kertoimen
	# string:in?
	# yks t�ytyy olla kokonaisluku, joka on
	# v�hint��n -1:n suuruinen. Pienemmill?
	# luvuilla tapahtuu jokin py�ristysvirhe.

	if (yks >= 0) {
		digit <- num %% 10 ^ (yks + 1)
		digit <- floor(digit / (10 ^ yks))
	} else {
		digit <- num * 10
		digit <- floor(digit %% 10)
	}
	digit <- as.character(digit)
	return(digit)
}