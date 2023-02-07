testFastaData <- function(inFile) {
  #  added by Lu Cheng, 11.11.2012
  if (!file.exists(inFile)) {
    stop('Fasta file ', inFile, ' does not exist!')
  }

  seqs <- load_fasta(inFile)
  heds <- colnames(seqs)
  ninds <- nrow(seqs) # not sure if this should be nrow, ncols or length

  data <- as.matrix(seqs)
  newData <- matrix(-9, nrow = nrow(data), ncol = ncol(data))
  newData[toupper(data) == 'A'] <- 1
  newData[toupper(data) == 'C'] <- 2
  newData[toupper(data) == 'G'] <- 3
  newData[toupper(data) == 'T'] <- 4

  data <- cbind(newData, seq(1, ninds))
  return(list("ninds" = ninds, "data" = data, "heds" = heds))
}
