#' @title Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.
#' @param data_matrix data_matrix
#' @param nclusters ncluster
#' @param Z Z

admixture_initialization <- function(data_matrix, nclusters, Z) {
  size_data <- size(data_matrix)
  nloci <- size_data[2] - 1
  n <- max(data_matrix[, ncol(data_matrix)])
  T <- cluster_own(Z, nclusters)
  initial_partition <- zeros(size_data[1], 1)
  for (i in 1:n) {
    kori <- T[i]
    here <- find(data_matrix[, ncol(data_matrix)] == i)
    for (j in 1:length(here)) {
      initial_partition[here[j], 1] <- kori
    }
  }
  return(initial_partition)
}
