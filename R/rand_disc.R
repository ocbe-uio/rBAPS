rand_disc <- function(CDF) {
	# %returns an index of a value from a discrete distribution using inversion method
	slump <- rand
	har <- find(CDF > slump)
	svar <- har(1)
	return(svar)
}
