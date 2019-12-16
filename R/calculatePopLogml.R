#' @title Calculate log marginal likelihood
#' @description Calculates fuzzy (log) marginal likelihood for a population of 
#' real values using estimate "fii" for the dispersion value, and Jeffreys prior
#' for the mean parameter.
#' @param points points
#' @param fii fii
calculatePopLogml <- function(points, fii) {
    n <- length(points)
    fuzzy_ones <- sum(points)
    fuzzy_zeros <- n - fuzzy_ones
    val = log(gamma(1)) -
        log(gamma(1 + n / fii)) +
        log(gamma(0.5 + fuzzy_ones / fii)) +
        log(gamma(0.5 + fuzzy_zeros / fii)) -
        log(gamma(0.5)) -
        log(gamma(0.5))
    return(val)
}

