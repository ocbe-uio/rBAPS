#' @title Clustering of individuals
#' @param data data file
#' @param format Format of the data c("FASTA", "VCF" ,"SAM", or "GenePop")
#' @importFrom utils read.delim
#' @export
greedyMix <- function(data, format) {
	format <- tolower(format)
	if (format == "fasta") {
		out <- load_fasta(data)
	} else if (format == "vcf") {
		stop("VCF files not yet supported." )
		# TODO #17: implement load_vcf()
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