newGetDistances <- function(data, rowsFromInd) {
	ninds <- max(data[, ncol(data)])
	nloci <- size(data, 2) - 1
	riviLkm <- choose(ninds, 2)

	empties <- find(data < 0)
	data[empties] <- 0
	data <- apply(data, 2, as.numeric) # max(noalle) oltava <256

	pariTaulu <- zeros(riviLkm, 2)
	aPointer <- 1
	for (a in 1:(ninds - 1)) {
		pariTaulu_rows <- aPointer:(aPointer + ninds - 1 - a)
		pariTaulu[pariTaulu_rows, 1] <- ones(ninds - a, 1) * a
		pariTaulu[pariTaulu_rows, 2] <- t((a + 1):ninds)
		aPointer <- aPointer + ninds - a
	}

	eka <- pariTaulu[, ones(1, rowsFromInd)]
	eka <- eka * rowsFromInd
	miinus <- repmat((rowsFromInd - 1):0, c(riviLkm, 1))
	eka <- eka - miinus

	toka <- pariTaulu[, ones(1, rowsFromInd) * 2]
	toka <- toka * rowsFromInd
	toka <- toka - miinus

	summa <- zeros(riviLkm, 1)
	vertailuja <- zeros(riviLkm, 1)

	rm(pariTaulu, miinus)

	x <- zeros(size(eka))
	x <- apply(x, 2, as.integer)
	y <- zeros(size(toka))
	y <- apply(y, 2, as.integer)

	for (j in 1:nloci) {
		for (k in 1:rowsFromInd) {
			x[, k] <- data[eka[, k], j]
			y[, k] <- data[toka[, k], j]
		}
		for (a in 1:rowsFromInd) {
			for (b in 1:rowsFromInd) {
				vertailutNyt <- as.double(x[, a] > 0 & y[, b] > 0)
				vertailuja <- vertailuja + vertailutNyt
				lisays <- (x[, a] != y[, b] & vertailutNyt)
				summa <- summa + as.double(lisays)
			}
		}
	}

	rm(x, y, vertailutNyt)
	nollat <- find(vertailuja == 0)
	dist <- zeros(length(vertailuja), 1)
	dist[nollat] <- 1
	muut <- find(vertailuja > 0)
	dist[muut] <- summa[muut] / vertailuja[muut]
	rm(summa, vertailuja)
	Z <- linkage(t(dist))
	return(list(Z = Z, dist = dist))
}