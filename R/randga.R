#' @title Generates random number from a Gamma distribution
#' @description Generates one random number from shape parameter a and rate parameter b
#' @param a shape
#' @param b rate
#' @return One realization of Gamma(a, b)
#' @details The generated random variable has mean a / b. It will be positively-skewed for small values, but converges to a symmetric distribution for very large numbers of a and b.
randga <- function(a, b) {
  flag <- 0
  if (a > 1) {
    c1 <- a - 1
    c2 <- (a - (1 / (6 * a))) / c1
    c3 <- 2 / c1
    c4 <- c3 + 2
    c5 <- 1 / sqrt(a)
    U1 <- 1
    while (flag == 0) {
      if (a <= 2.5) {
        U1 <- rand()
        U2 <- rand()
      } else {
        while (!(U1 > 0 & U1 < 1)) {
          U1 <- rand()
          U2 <- rand()
          U1 <- U2 + c5 * (1 - 1.86 * U1)
        }
      }
      W <- c2 * U2 / U1
      if (c3 * U1 + W + (1 / W) <= c4) {
        flag <- 1
        g <- c1 * W / b
      } else if (c3 * log(U1) - log(W) + W < 1) {
        flag <- 1
        g <- c1 * W / b
      } else {
        U1 <- -1
      }
    }
  } else if (a == 1) {
    g <- sum(-(1 / b) * log(rand(a, 1)))
  } else {
    while (flag == 0) {
      U <- rand(2, 1)
      if (U[1] > exp(1) / (a + exp(1))) {
        g <- -log(((a + exp(1)) * (1 - U[1])) / (a * exp(1)))
        if (U[2] <= g^(a - 1)) {
          flag <- 1
        }
      } else {
        g <- ((a + exp(1)) * U[1] / ((exp(1))^(1 / a)))
        if (U[2] <= exp(-g)) {
          flag <- 1
        }
      }
    }
    g <- g / b
  }
  return(g)
}
