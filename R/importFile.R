#' @title Import data file
#' @description Imports data from several formats (FASTA, VCF, SAM, BAM,
#' Genepop).
#' @param data raw dataset
#' @param format data format (guesses from extension if not provided)
#' @param verbose if \code{TRUE}, prints extra output information
#' @return The data in a format that can be used by the other functions
#' @export
#' @examples
#' path_inst <- system.file("extdata", "", package = "rBAPS")
#' importFile(file.path(path_inst, "FASTA_clustering_haploid.fasta"))
importFile <- function(data, format, verbose) {
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
