baps.globals <- new.env(parent = emptyenv())

assign("COUNTS", array(0, dim = c(0, 0, 0)), envir = baps.globals)
assign("SUMCOUNTS", array(0, dim = c(0, 0)), envir = baps.globals)
assign("PARTITION", array(1, dim = 0), envir = baps.globals)
assign("POP_LOGML", array(1, dim = 0), envir = baps.globals)
assign("LOGDIFF", array(1, dim = c(0, 0)), envir = baps.globals)
# If handling globas break, try other ideas from
# https://stackoverflow.com/a/65252740/1169233 and
# https://stackoverflow.com/questions/12598242/
