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
#' greedyMix(data)
greedyMix <- function(data, format, verbose = TRUE) {
  # Parsing data format ------------------------------------------------------

  if (missing(format)) {
    format <- gsub(".*\\.(.+)$", "\\1", data)
    message("Format not provided. Guessing from file extension: ", format)
  }
  format <- tolower(format)

  # Dispatching to proper loading function -----------------------------------

  if (format == "fasta") {
    out <- load_fasta(data)
  } else if (format == "vcf") {
    out <- vcfR::read.vcfR(data, verbose = verbose)
  } else if (format == "sam") {
    stop(
      "SAM files not directly supported. ",
      "Install the samtools software and execute\n\n",
      "samtools view -b ", data, " > out_file.bam\n\nto convert to BAM ",
      "and try running this function again with 'format=BAM'"
    )
  } else if (format == "bam") {
    out <- Rsamtools::scanBam(data)
  } else if (format == "genepop") {
    if (toupper(adegenet::.readExt(data)) == "TXT") {
      message("Creating a copy of the file with the .gen extension")
      dataGen <- gsub("txt", "gen", data)
      file.copy(data, dataGen)
      out <- adegenet::read.genepop(dataGen)
    } else {
      out <- adegenet::read.genepop(data)
    }
  } else {
    stop("Format not supported.")
  }
  return(out)
}
