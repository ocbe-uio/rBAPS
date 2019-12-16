#' @title Learn partition (modified)
#' @param ordered ordered
#' @return part
#' @description This function is called only if some individual has less than 
#' 90 per cent non-missing data. The function uses fuzzy clustering for the
#' "non-missingness" values, finding maximum three clusters. If two of the
#' found clusters are such that all the values are >0.9, then those two are 
#' further combined.
learn_partition_modified <- function(ordered) {
    part <- learn_simple_partition(ordered, 0.05)
    # TODO: find learn_simple_partition
    nclust <- length(unique(part))
    if (nclust == 3) {
        mini_1 <- min(ordered(which(part == 1))) # ASK: what is ordered()?
        mini_2 <- min(ordered(which(part == 2)))
        mini_3 <- min(ordered(which(part == 3)))
        if (mini_1 > 0.9 & mini_2 > 0.9) {
            part[part == 2] <- 1
            part[part == 3] <- 2
        } else if (mini_1 > 0.9 & mini_3 > 0.9) {
            part[part == 3] <- 1
        } else if (mini_2 > 0.9 & mini_3 > 0.9) {
            # This is the one happening in practice, since the values are
            # ordered, leading to mini_1 <= mini_2 <= mini_3
            part[part == 3] <- 2
        }
    }
    return(part)
}