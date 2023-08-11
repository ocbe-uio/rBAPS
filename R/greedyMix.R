#' @title Clustering of individuals
#' @param data data file
#' @param format Data format. Format supported: "FASTA", "VCF" ,"BAM", "GenePop"
#' @param partitionCompare a list of partitions to compare
#' @param ninds number of individuals
#' @param rowsFromInd a list of rows for each individual
#' @param noalle number of alleles
#' @param adjprior ajuster prior probabilities
#' @param npops number of populations
#' @param priorTerm prior terms
#' @param counts counts
#' @param sumcounts sumcounts
#' @param max_iter maximum number of iterations
#' @param alleleCodes allele codes
#' @param inp input file
#' @param popnames population names
#' @param fixedK if \code{TRUE}, the number of populations is fixed
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
  data, format, partitionCompare = NULL, ninds = 1L, npops = 1L,
  counts = NULL, sumcounts = NULL, max_iter = 100L, alleleCodes = NULL,
  inp = NULL, popnames = NULL, fixedK = FALSE, verbose = FALSE
) {
  # Importing and handling data ================================================
  data <- importFile(data, format, verbose)
  data <- handleData(data, tolower(format))
  c <- list(
    noalle = data[["noalle"]],
    data = data[["newData"]],
    adjprior = data[["adjprior"]],
    priorTerm = data[["priorTerm"]],
    rowsFromInd = data[["rowsFromInd"]]
  )

  # Comparing partitions =======================================================
  if (!is.null(partitionCompare)) {
    logmls <- comparePartitions(
      c[["data"]], nrow(c[["data"]]), partitionCompare[["partitions"]], ninds,
      c[["rowsFromInd"]], c[["noalle"]], c[["adjprior"]]
    )
  }


  # Generating partition summary ===============================================
  ekat <- seq(1L, c[["rowsFromInd"]], ninds * c[["rowsFromInd"]]) # ekat = (1:rowsFromInd:ninds*rowsFromInd)';
  c[["rows"]] <- c(ekat, ekat + c[["rowsFromInd"]] - 1L) # c.rows = [ekat ekat+rowsFromInd-1]
  logml_npops_partitionSummary <- indMixWrapper(c, npops, counts, sumcounts, max_iter, fixedK, verbose);
  logml <- logml_npops_partitionSummary[["logml"]]
  npops <- logml_npops_partitionSummary[["npops"]]
  partitionSummary <- logml_npops_partitionSummary[["partitionSummary"]]

  # Generating output object ===================================================
  out <- list(
    "alleleCodes" = alleleCodes, "adjprior" = c[["adjprior"]],
    "popnames" = popnames, "rowsFromInd" = c[["rowsFromInd"]],
    "data" = c[["data"]], "npops" = npops, "noalle" = c[["noalle"]],
    "mixtureType" = "mix", "logml" = logml
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
