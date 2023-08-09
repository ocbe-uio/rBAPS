#' @title Clustering of individuals
#' @param data data file
#' @param format Data format. Format supported: "FASTA", "VCF" ,"BAM", "GenePop"
#' @param verbose if \code{TRUE}, prints extra output information
#' @importFrom utils read.delim
#' @importFrom vcfR read.vcfR
#' @importFrom Rsamtools scanBam
#' @importFrom adegenet read.genepop .readExt
#' @references Samtools: a suite of programs for interacting
#' with high-throughput sequencing data. <http://www.htslib.org/>
#' @export
#' @examples
#' data <- system.file("extdata", "FASTA_clustering_haploid.fasta", package = "rBAPS")
#' greedyMix(data, "fasta")
greedyMix <- function(
  data, format, partitionCompare, ninds, rowsFromInd, noalle,
  adjprior, priorTerm, alleleCodesinp, popnames, fixedK = FALSE,
  partition_compare = FALSE, verbose = TRUE
) {
  # Importing and handling data ================================================
  data <- importFile(data, format, verbose)

  if (partition_compare) {
    logmls <- comparePartitions(
      data, nrows(data), partitionCompare[["partitions"]], ninds, rowsFromInd,
      noalle, adjprior
    )
  }


  # Generating partition summary ===============================================
  logml_npops_partitionSummary <- indMixWrapper(data, npops, fixedK);
  logml <- logml_npops_partitionSummary[["logml"]]
  npops <- logml_npops_partitionSummary[["npops"]]
  partitionSummary <- logml_npops_partitionSummary[["partitionSummary"]]
  stopifnot(logml != 1)

  # Writing mixture info =======================================================
  changesInLogml <- writeMixtureInfo(
    logml, rowsFromInd, data, adjprior, priorTerm, NULL, inp, partitionSummary,
    popnames, fixedK
  )

  # Returning results ==========================================================
  return(
    list(
      "alleleCodes" = alleleCodes, "adjprior" = adjprior, "popnames" = popnames,
      "rowsFromInd" = rowsFromInd, "data" = data, "npops" = npops,
      "noalle" = noalle, "mixtureType" = "mix", "logml" = logml,
      "changesInLogml" = changesInLogml
    )
  )

}
