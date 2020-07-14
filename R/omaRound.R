omaRound <- function(num) {
	# Py�rist�� luvun num 1 desimaalin tarkkuuteen
	num <- num * 10
	num <- round(num)
	num2 <- num / 10
	return(num2)
}