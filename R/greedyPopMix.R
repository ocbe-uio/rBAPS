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
#' @references Samtools: a suite of programs for interacting
#' with high-throughput sequencing data. <http://www.htslib.org/>
#' @export
greedyPopMix <- function(data, format, partitionCompare = NULL, verbose = TRUE) {
  # Replacing original file reading code with greedyMix()
  greedyMix(data, format, verbose)
  # TODO: find out where the elements above come from. Maybe greedyMix should create them?
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
        partitionInd <- zeros(rows(end), 1)
        partitionSample <- partitions[, i]
        for (j in 1:nsamplingunits) {
          partitionInd[c$rows[j, 1]:c.rows[j, 2]] <- partitionSample[j]
        }
        partitionLogml[i] <- initialCounts(
          partitionInd, data[, 1:(ncol(data) - 1)], npops, c$rows, noalle,
          adjprior
        )
    }
    # return the logml result
    partitionCompare$logmls <- partitionLogml
  }
  data = data(:,1:end-1);
  data <- data[, 1:(ncol(data) - 1)]
  changesInLogml <- writeMixtureInfoPop(
    logml, rows, data, adjprior, priorTerm,
    outp, inp, partitionSummary, popnames, fixedK
  )
  talle <- questdlg(
    'Do you want to save the mixture populations so that you can use them later in admixture analysis?',
    'Save results?', c('Yes', 'No'), 'Yes'
  )
  if (isequal(talle, 'Yes')) {
    waitALittle()
    filename_pathname <- uiputfile()
    if (rowsFromInd == 0) {
      # BAPS format was used, rowsFromInd is not known.
      popnames_rowsFromInd <- findOutRowsFromInd(popnames, rows)
      popnames_rowsFromInd$popnames -> popnames
      popnames_rowsFromInd$rows -> rows
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
    c$mixtureType = 'popMix'
    c$groupPartition <- groupPartition
    c$rows <- rows
    c$logml <- logml
    c$changesInLogml <- changesInLogml
  }
}
