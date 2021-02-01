COUNTS <- array(0, dim=c(100, 100, 100))
SUMCOUNTS <- array(0, dim=c(100, 100))
PARTITION <- vector()
POP_LOGML <- vector()
LOGDIFF <- vector()
# If handling globas break, try other ideas from https://stackoverflow.com/a/65252740/1169233


#' @import utils
utils::globalVariables(
	c("PARTITION", "COUNTS", "SUMCOUNTS", "LOGDIFF", "POP_LOGML", "GAMMA_LN")
)