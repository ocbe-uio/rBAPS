#' load_fasta
#'
#' Loads a fasta file into matrix format ready for
#' running the hierBAPS algorithm.
#'
#' @param msa Either the location of a fasta file or ape DNAbin object containing the multiple sequence alignment data to be clustered
#' @param keep_singletons A logical indicating whether to consider singleton mutations in calculating the clusters
#' @param output_numbers A logical indicating whether to output the data as
#' numbers (TRUE) or letters (FALSE)
#'
#' @return A character matrix with filtered SNP data
#'
#' @examples
#' msa <- system.file("extdata", "seqs.fa", package = "rBAPS")
#' snp.matrix <- rBAPS:::load_fasta(msa)
#' @author Gerry Tonkin-Hill, Waldir Leoncio
#' @seealso rhierbaps::load_fasta
#' @importFrom ape read.FASTA as.DNAbin
load_fasta <- function(msa, keep_singletons = FALSE, output_numbers = TRUE) {

  # Check inputs
  if (is(msa, "character")) {
    if (!file.exists(msa)) stop("Invalid msa or the file does not exist!")
    seqs <- ape::read.FASTA(msa)
  } else if (is(msa, "matrix")) {
    seqs <- ape::as.DNAbin(msa)
  } else if (is(msa, "DNAbin")) {
    seqs <- msa
  } else {
    stop("incorrect input for msa!")
  }
  if (!is.logical(keep_singletons)) {
    stop("Invalid keep_singletons! Must be one of TRUE/FALSE.")
  }

  # Load sequences using ape. This does a lot of the checking for us.
  seq_names <- labels(seqs)
  seqs <- as.character(as.matrix(seqs))
  rownames(seqs) <- seq_names
  seqs[is.na(seqs)] <- "-"

  # Validation -----------------------------------------------------------------
  if (nrow(seqs) < 3) stop("Less than 3 sequences!")
  if (any(!(as.vector(tolower(seqs)) %in% c("a", "c", "g", "t", "n", "-")))) {
    warning("Characters not in acgtnACGTN- will be treated as missing (-)...")
  }

  # Remove conserved columns
  conserved <- colSums(t(t(seqs) == seqs[1, ])) == nrow(seqs)
  seqs <- seqs[, !conserved]

  if (!keep_singletons) {
    # remove_singletons as they are uninformative in the algorithm
    is_singleton <- apply(seqs, 2, function(x) {
      tab <- table(x)
      return(x %in% names(tab)[tab == 1])
    })
    seqs[is_singleton] <- "-"
  }

  # Convert gaps and unknowns to same symbol
  seqs[seqs == "n"] <- "-"

  # Replace letters with numbers, dashes with zeros
  if (output_numbers) {
    seqs <- matrix(match(seqs, c("a", "c", "g", "t")), nrow(seqs))
    seqs[is.na(seqs)] <- 0
  }

  return(seqs)
}
