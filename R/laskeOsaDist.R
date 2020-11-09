laskeOsaDist <- function(inds2, dist, ninds) {
	# % Muodostaa dist vektorista osavektorin, joka sis�lt�� yksil�iden inds2
	# % v�liset et�isyydet. ninds=kaikkien yksil�iden lukum��r�.

	ninds2 <- length(inds2)
	apu <- zeros(nchoosek(ninds2, 2), 2)
	rivi <- 1
	for (i in 1:ninds2-1) {
		for (j in i+1:ninds2) {
			apu[rivi, 1] <- inds2[i]
			apu[rivi, 2] <- inds2[j]
			rivi <- rivi + 1
		}
	}
	apu <- (apu[, 1]-1) * ninds - apu[, 1] / 2 *
		(apu[, 1]-1) + (apu[, 2] - apu[, 1])
	dist2 <- dist(apu)
	return(dist2)
}