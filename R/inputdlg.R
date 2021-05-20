#' @title Gather user input
#' @description Replicates the functionality of the homonymous function in Matlab (sans dialog box)
#' @param prompt Text field with user instructions
#' @param dims number of dimensions in the answwers
#' @param definput default value of the input
#' @export
inputdlg <- function(prompt, dims=1, definput=NULL) {
    if (!is.null(definput)) {
        prompt <- append(prompt, paste0(" (default: ", definput, ")"))
    }
    input_chr <- readline(paste0(prompt, ": "))
    if (input_chr == "") input_chr <- definput
    input_chr_or_num <- tryCatch(
        as.numeric(input_chr), warning = function(w) input_chr
    )
    return(input_chr_or_num)
}