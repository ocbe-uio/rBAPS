#' @title Clustering of individuals
#' @param data data file
#' @param format Format of the data c("FASTA", "VCF" ,"SAM", or "GenePop")
#' @param verbose if \code{TRUE}, prints extra output information
#' @importFrom utils read.delim
#' @export
greedyMix <- function(data, format, verbose = TRUE) {
	format <- tolower(format)
	if (format == "fasta") {
		out <- load_fasta(data)
	} else if (format == "vcf") {
		# TODO #17: implement load_vcf()
		out <- vcfR::read.vcfR(data, verbose = verbose)
	} else if (format == "sam") {
		stop("SAM files not yet supported." )
		# TODO #18: implement load_sam()
	} else if(format == "genepop") {
		# TODO #19: implement load_genepop()
		stop("GenePop files not yet supported." )
	} else {
		stop("Format not supported.")
	}
	return(out)
}