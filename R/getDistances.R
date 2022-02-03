getDistances <- function(data_matrix, nclusters) {

  # %finds initial admixture clustering solution with nclusters clusters, uses simple mean Hamming distance
  # %gives partition in 8 - bit format
  # %allocates all alleles of a single individual into the same basket
  # %data_matrix contains #Loci + 1 columns, last column indicate whose alleles are placed in each row,
  # %i.e. ranges from 1 to #individuals. For diploids there are 2 rows per individual, for haploids only a single row
  # %missing values are indicated by zeros in the partition and by negative integers in the data_matrix.

  size_data <- size(data_matrix)
  nloci <- size_data[2] - 1
  n <- base::max(data_matrix[, ncol(data_matrix)])
  distances <- zeros(choose(n, 2), 1)
  pointer <- 1
  for (i in 1:n - 1) {
    i_data <- data_matrix[
      matlab2r::find(data_matrix[, ncol(data_matrix)] == i),
      1:nloci
    ]
    for (j in (i + 1):n) {
      d_ij <- 0
      j_data <- data_matrix[matlab2r::find(data_matrix[, ncol()] == j), 1:nloci]
      vertailuja <- 0
      for (k in 1:size(i_data, 1)) {
        for (l in 1:size(j_data, 1)) {
          here_i <- matlab2r::find(i_data[k, ] >= 0)
          here_j <- matlab2r::find(j_data[l, ] >= 0)
          here_joint <- intersect(here_i, here_j)
          vertailuja <- vertailuja + length(here_joint)
          d_ij <- d_ij + length(
            matlab2r::find(i_data[k, here_joint] != j_data[l, here_joint])
          )
        }
      }
      d_ij <- d_ij / vertailuja
      distances[pointer] <- d_ij
      pointer <- pointer + 1
    }
  }

  Z <- linkage(t(distances))
  return(list(Z = Z, distances = distances))
}
