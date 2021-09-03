#' @title Clustering of individuals
#' @param data data file
#' @param format Format of the data c("FASTA", "VCF" ,"BAM", or "GenePop")
#' @param verbose if \code{TRUE}, prints extra output information
#' @importFrom utils read.delim
#' @importFrom vcfR read.vcfR
#' @importFrom Rsamtools scanBam
#' @references Samtools: a suite of programs for interacting
#' with high-throughput sequencing data. <http://www.htslib.org/>
#' @export
greedyMix <- function(data, format, verbose = TRUE) {
	format <- tolower(format)
	if (format == "fasta") {
		out <- load_fasta(data)
	} else if (format == "vcf") {
		out <- vcfR::read.vcfR(data, verbose = verbose)
	} else if (format == "sam") {
		stop(
			"SAM files not directly supported. ",
			"Install the samtools software and execute ",
			"'samtools view -b in_file.sam > out_file.bam' to convert to BAM ",
			"and try running this function again with 'format=BAM'"
		)
	} else if (format == "bam") {
		out <- Rsamtools::scanBam(data)
	} else if (format == "genepop") {
		# TODO #19: implement load_genepop()
		stop("GenePop files not yet supported." )
	} else {
		stop("Format not supported.")
	}
	return(out)
}