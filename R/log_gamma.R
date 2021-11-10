#' @title Log Gamma
#' @description Equal to log(gamma(x)) with special handling of x < 0 for
#' Matlab compatibility
#' @param x number
#' @return log(gamma(x)) for x > 0, Inf otherwise
log_gamma <- function(x) {
  ifelse(x > 0, log(gamma(x)), Inf)
}
