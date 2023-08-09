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
  data, format, partitionCompare = NULL, ninds = NULL, rowsFromInd = NULL,
  noalle = NULL, adjprior = NULL, npops = 1L, priorTerm = NULL,
  alleleCodes = NULL, inp = NULL, popnames = NULL, fixedK = FALSE, verbose = TRUE
) {
  # Importing and handling data ================================================
  data <- importFile(data, format, verbose)
  c <- list(
    # TODO: get elements from handleData()?
    noalle = noalle,
    rows = NA,
    data = data,
    adjprior = adjprior,
    priorTerm = priorTerm,
    rowsFromInd = rowsFromInd
  )

  # Comparing partitions =======================================================
  if (!is.null(partitionCompare)) {
    logmls <- comparePartitions(
      data, nrows(data), partitionCompare[["partitions"]], ninds, rowsFromInd,
      noalle, adjprior
    )
  }


  # Generating partition summary ===============================================
  logml_npops_partitionSummary <- indMixWrapper(c, npops, fixedK, verbose);
  logml <- logml_npops_partitionSummary[["logml"]]
  npops <- logml_npops_partitionSummary[["npops"]]
  partitionSummary <- logml_npops_partitionSummary[["partitionSummary"]]

  # Generating output object ===================================================
  out <- list(
      "alleleCodes" = alleleCodes, "adjprior" = adjprior, "popnames" = popnames,
      "rowsFromInd" = rowsFromInd, "data" = data, "npops" = npops,
      "noalle" = noalle, "mixtureType" = "mix", "logml" = logml
    )
  if (logml == 1) {
    return(out)
  }

  # Writing mixture info =======================================================
  changesInLogml <- writeMixtureInfo(
    logml, rowsFromInd, data, adjprior, priorTerm, NULL, inp, partitionSummary,
    popnames, fixedK
  )

  # Updateing results ==========================================================
  return(c(out, "changesInLogml" = changesInLogml))
}
