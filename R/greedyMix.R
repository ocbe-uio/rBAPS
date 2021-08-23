#' @title Clustering of individuals
#' @param tietue File
#' @param format Format of the data ("BAPS", "GenePop" or "Preprocessed")
#' @param savePreProcessed Save the pre-processed data?
#' @param filePreProcessed Is the file already processed?
#' @importFrom utils read.delim
#' @export
greedyMix <- function(
	tietue,
	format           = NULL,
	savePreProcessed = NULL,
	filePreProcessed = NULL
) {
	stop(
		"greedyMix() has been superseded by load_fasta().",
		" Please change your code to use the latter instead of the former.",
		" If you believe the error is internal to rBAPS, please open",
		" a new issue (link on the package DESCRIPTION file)."
		)
}
