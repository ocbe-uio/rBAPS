#' @title Convert from FASTA to BAPS
#' @description Converts a file (not an R object) from FASTA to BAPS format
#' @param file filename of FASTA file
#' @return `data` in BAPS format
#' @author Waldir Leoncio
#' @export
#' @examples
#' file <- system.file("extdata", "FASTA_clustering_haploid.fasta", package = "rBAPS")
#' convert_FASTA_to_BAPS(file)
convert_FASTA_to_BAPS <- function(file) {
  data <- load_fasta(file) # Processing data
  data <- cbind(data, seq_len(nrow(data))) # Add IDs of individuals (sequential)
  data[data == 0] <- -9 # Because zeros (missing) in BAPS are coded as -9
  return(data)
}
