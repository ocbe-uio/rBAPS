COUNTS <- array(0, dim=c(100, 100, 100))
SUMCOUNTS <- array(0, dim=c(100, 100))
PARTITION <- array(1, dim=c(100))
POP_LOGML <-  array(1, dim=c(100))
LOGDIFF <-  array(1, dim=c(100, 100))
# If handling globas break, try other ideas from https://stackoverflow.com/a/65252740/1169233


#' @import utils
utils::globalVariables(
	c("PARTITION", "COUNTS", "SUMCOUNTS", "LOGDIFF", "POP_LOGML", "GAMMA_LN")
)