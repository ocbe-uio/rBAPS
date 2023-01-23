testFastaData <- function(inFile) {
  #  added by Lu Cheng, 11.11.2012
  if (!exists(inFile, 'file')) {
    stop('Fasta file ', inFile, ' does not exist!')
  }

  seqs <- load_fasta(inFile)
  heds <- colnames(seqs)
  ninds <- length(seqs)

  data <- as.matrix(seqs)
  newData <- ones(size(data)) * -9
  newData[toupper(data) == 'A'] <- 1
  newData[toupper(data) == 'C'] <- 2
  newData[toupper(data) == 'G'] <- 3
  newData[toupper(data) == 'T'] <- 4
  data <- c(newData, t(1:ninds))
  return(list("ninds" = ninds, "data" = data, "heds" = heds))
}
