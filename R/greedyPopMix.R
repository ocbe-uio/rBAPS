#' @title Clustering of pop individuals
#' @param data data file
#' @param format Data format. Format supported: "FASTA", "VCF" ,"BAM", "GenePop"
#' @param verbose if \code{TRUE}, prints extra output information
#' @param partitionCompare a properly-named list. Proper names include
#' "partitions"
#' @importFrom utils read.delim
#' @importFrom vcfR read.vcfR
#' @importFrom Rsamtools scanBam
#' @importFrom adegenet read.genepop .readExt
#' @importFrom matlab2r uiputfile
#' @references Samtools: a suite of programs for interacting
#' with high-throughput sequencing data. <http://www.htslib.org/>
greedyPopMix <- function(data, format, partitionCompare = NULL, verbose = TRUE
) {
  # Replacing original file reading code with greedyMix()
  rawdata <- greedyMix(data, format, verbose)

  # Other function calls to produce necessary objects
  data_greedyMix_handle <- handlePopData(rawdata)
  data <- data_greedyMix_handle$data
  rowsFromInd <- data_greedyMix_handle$rowsFromInd
  alleleCodes <- data_greedyMix_handle$alleleCodes
  noalle <- data_greedyMix_handle$noalle
  adjprior <- data_greedyMix_handle$adjprior
  priorTerm <- data_greedyMix_handle$priorTerm
  rm(data_greedyMix_handle)
  Z_dist <- getPopDistancesByKL(adjprior)
  Z <- Z_dist$Z
  dist <- Z_dist$dist
  rm(Z_dist)
  a_data <- data[, 1:(ncol(data) - 1)]
  sumcounts_counts_logml <- initialPopCounts(a_data, npops, rows, noalle, adjprior)
  logml <- sumcounts_counts_logml$logml
  rm(sumcounts_counts_logml)
  c <- list()
  c$data <- data
  c$rows <- rows
  c$alleleCodes <- alleleCodes
  c$noalle <- noalle
  c$adjprior <- adjprior
  c$priorTerm <- priorTerm
  c$dist <- dist
  c$Z <- Z
  c$rowsFromInd <- rowsFromInd
  #  partition compare
  if (!is.null(partitionCompare)) {
    nsamplingunits <- size(rows, 1)
    partitions <- partitionCompare$partitions
    npartitions <- size(partitions, 2)
    partitionLogml <- zeros(1, npartitions)
    for (i in 1:npartitions) {
      # number of unique partition lables
        npops <- length(unique(partitions[, i]))
        partitionInd <- zeros(length(rows), 1)
        partitionSample <- partitions[, i]
        for (j in 1:nsamplingunits) {
          partitionInd[c$rows[j, 1]:c$rows[j, 2]] <- partitionSample[j]
        }
        partitionLogml[i] <- initialCounts(
          partitionInd, data[, 1:(ncol(data) - 1)], npops, c$rows, noalle,
          adjprior
        )
    }
    # return the logml result
    partitionCompare$logmls <- partitionLogml
  }
  data <- data[, 1:(ncol(data) - 1)]
  logml_npops_partitionSummary <- indMix(c)
  logml_npops_partitionSummary$logml -> logml
  logml_npops_partitionSummary$npops -> npops
  logml_npops_partitionSummary$partitionSummary -> partitionSummary
  rm(logml_npops_partitionSummary)
  changesInLogml <- writeMixtureInfoPop(
    logml, rows, data, adjprior, priorTerm,
    NULL, NULL, partitionSummary, popnames, fixedK = FALSE
  )
  talle <- questdlg(
    "Do you want to save the mixture populations so that you can use them later in admixture analysis?",
    "Save results?", c("Yes", "No"), "Yes"
  )
  if (tolower(talle) == "yes") {
    waitALittle()
    filename_pathname <- uiputfile()
    if (rowsFromInd == 0) {
      # BAPS format was used, rowsFromInd is not known.
      popnames_rowsFromInd <- findOutRowsFromInd(popnames, rows)
      popnames <- popnames_rowsFromInd$popnames
      rows <- popnames_rowsFromInd$rows
      rm(popnames_rowsFromInd)
    }
    groupPartition <- PARTITION
    fiksaaPartitioYksiloTasolle(rows, rowsFromInd)
    c$PARTITION <- PARTITION
    c$COUNTS <- COUNTS
    c$SUMCOUNTS <- SUMCOUNTS
    c$alleleCodes <- alleleCodes
    c$adjprior <- adjprior
    c$rowsFromInd <- rowsFromInd
    c$popnames <- popnames
    c$data <- data
    c$npops <- npops
    c$noalle <- noalle
    c$mixtureType <- "popMix"
    c$groupPartition <- groupPartition
    c$rows <- rows
    c$logml <- logml
    c$changesInLogml <- changesInLogml
  }
}
