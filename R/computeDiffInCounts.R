computeDiffInCounts <- function(rows, max_noalle, nloci, data) {
	# % Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
	# % lukum��r�t (vastaavasti kuin COUNTS:issa), jotka ovat data:n
	# % riveill� rows. rows pit�� olla vaakavektori.

	diffInCounts <- zeros(max_noalle, nloci)
	for (i in seq_len(nrow(data))) {
		row <- data[i, ]
		notEmpty <- as.matrix(find(row>=0))

		if (length(notEmpty) > 0) {
			diffInCounts[row(notEmpty) + (notEmpty - 1) * max_noalle] <-
				diffInCounts[row(notEmpty) + (notEmpty - 1) * max_noalle] + 1
		}
	}
	diffInCounts <- matrix(
		data = diffInCounts[!is.na(diffInCounts)],
		nrow = max_noalle,
		ncol = nloci,
		byrow = TRUE
	)
	return(diffInCounts)
}
