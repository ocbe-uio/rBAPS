findOutRowsFromInd <- function(popnames, rows, ploidisuus = NULL) {
  if (is.null(ploidisuus)) {
    ploidisuus <- questdlg(
      quest = 'Specify the type of individuals in the data',
      dlgtitle = 'Individual type?',
      btn = c('Haploid', 'Diploid', 'Tetraploid'),
      defbtn = 'Diploid'
    )
  }

  rowsFromInd <- switch(ploidisuus,
    'Haploid' = 1,
    'Diploid' = 2,
    'Tetraploid' = 4
  )

  popnames2 <- popnames * NA
  if (!is.null(popnames)) {
      for (i in seq_len(size(rows, 1))) {
          popnames2[i, 1] <- popnames[i, 1]
          rivi <- rows[i, 1]:rows[i, 2]
          popnames2[i, 2] <- rivi[rowsFromInd] / rowsFromInd
      }
  }
	return(list(popnames2 = popnames2, rowsFromInd = rowsFromInd))
}
