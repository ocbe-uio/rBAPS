computeDiffInCounts <- function(rows, max_noalle, nloci, data) {
  # % Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
  # % lukum��r�t (vastaavasti kuin COUNTS:issa), jotka ovat data:n
  # % riveill� rows. rows pit�� olla vaakavektori.

  diffInCounts <- zeros(max_noalle, nloci)
  for (i in rows) { # yep, just one iteration
    row <- data[i, ]
    notEmpty <- as.matrix(matlab2r::find(row >= 0))

    if (length(notEmpty) > 0) {
      element <- row[notEmpty] + (notEmpty - 1) * max_noalle
      diffInCounts[element] <- diffInCounts[element] + 1
    }
  }
  return(diffInCounts)
}
