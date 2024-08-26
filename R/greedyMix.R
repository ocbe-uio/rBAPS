#' @title Clustering of individuals
#' @param data data file
#' @param format Data format. Format supported: "FASTA", "VCF" ,"BAM", "GenePop"
#' @param partitionCompare a list of partitions to compare
#' @param npops number of populations
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
#' data <- system.file("extdata", "BAPS_clustering_diploid.txt", package = "rBAPS")
#' greedyMix(data, "baps")
greedyMix <- function(
  data, format = gsub("^.*\\.", "", data), partitionCompare = NULL, npops = 3L,
  counts = NULL, sumcounts = NULL, max_iter = 100L, alleleCodes = NULL,
  inp = NULL, popnames = NULL, fixedK = FALSE, verbose = FALSE
) {
  # Importing and handling data ================================================
  if (tolower(format) %in% "fasta") {
    data <- convert_FASTA_to_BAPS(data)
    format <- "baps"
  }
  if (tolower(format) %in% "baps") {
    data <- process_BAPS_data(data, NULL)
    c <- list(
      noalle = data[["noalle"]],
      data = data[["data"]],
      adjprior = data[["adjprior"]],
      priorTerm = data[["priorTerm"]],
      rowsFromInd = data[["rowsFromInd"]],
      Z = data[["Z"]],
      dist = data[["dist"]]
    )
  } else {
    data <- importFile(data, format, verbose)
    data <- handleData(data, tolower(format))
    c <- list(
      noalle = data[["noalle"]],
      data = data[["newData"]],
      adjprior = data[["adjprior"]],
      priorTerm = data[["priorTerm"]],
      rowsFromInd = data[["rowsFromInd"]],
      Z = data[["Z"]],
      dist = data[["dist"]]
    )
  }

  # Comparing partitions =======================================================
  ninds <- length(unique(c[["data"]][, ncol(c[["data"]])]))
  if (!is.null(partitionCompare)) {
    logmls <- comparePartitions(
      c[["data"]], nrow(c[["data"]]), partitionCompare[["partitions"]], ninds,
      c[["rowsFromInd"]], c[["noalle"]], c[["adjprior"]]
    )
  }

  # Generating partition summary ===============================================
  ekat <- seq(1L, ninds * c[["rowsFromInd"]], c[["rowsFromInd"]])
  c[["rows"]] <- cbind(ekat, ekat + c[["rowsFromInd"]] - 1L)
  logml_npops_partitionSummary <- indMixWrapper(c, npops, counts, sumcounts, max_iter, fixedK, verbose) # FIXME: not working for FASTA data
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
    logml, c[["rowsFromInd"]], c[["data"]], c[["adjprior"]], c[["priorTerm"]],
    NULL, inp, partitionSummary, popnames, fixedK, verbose
  )

  # Updateing results ==========================================================
  return(c(out, list("changesInLogml" = changesInLogml)))
}
