globals <- new.env(parent = emptyenv())

assign("COUNTS", array(0, dim = c(0, 0, 0)), envir = globals)
assign("SUMCOUNTS", array(0, dim = c(0, 0)), envir = globals)
assign("PARTITION", array(1, dim = 0), envir = globals)
assign("POP_LOGML", array(1, dim = 0), envir = globals)
assign("LOGDIFF", array(1, dim = c(0, 0)), envir = globals)
assign("CLOQCOUNTS", array(0, dim = c(0, 0)), envir = globals)
assign("SEPCOUNTS", array(0, dim = c(0, 0)), envir = globals)
assign("GAMMA_LN", array(0, dim = c(0, 0)), envir = globals)
# If handling globas break, try other ideas from
# https://stackoverflow.com/a/65252740/1169233 and
# https://stackoverflow.com/questions/12598242/
