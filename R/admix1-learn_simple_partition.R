#' @title Learn simple partition
#' @param ordered_points ordered_points
#' @param fii fii
#' @description Goes through all the ways to divide the points into two or 
#' three groups. Chooses the partition which obtains highest logml.
#' @export
learn_simple_partition <- function(ordered_points, fii) {
    npoints <- length(ordered_points)
    
    # One cluster:
    val <- calculatePopLogml(ordered_points, fii)
    bestValue <- val
    best_type <- 'single'
    
    # Two clusters:
    for (i in 1:(npoints - 1)) {
        # The right endpoint of the first cluster.
        val_1 <- calculatePopLogml(ordered_points[1:i], fii)
        val_2 <- calculatePopLogml(ordered_points[(i + 1):length(ordered_points)], fii)
        total <- val_1 + val_2
        if (total > bestValue) {
            bestValue <- total
            best_type <- 'double'
            best_i <- i
        }
    }

    # Three clusters:
    for (i in 1:(npoints - 2)) {
        for (j in (i + 1):(npoints - 1)) {
            val_1 <- calculatePopLogml(ordered_points[1:i], fii)
            val_2 <- calculatePopLogml(ordered_points[(i + 1):j], fii)
            val_3 <- calculatePopLogml(ordered_points[(j + 1):length(ordered_points)], fii)
            total <- val_1 + val_2 + val_3
            if (total > bestValue) {
                bestValue <- total
                best_type <- 'triple'
                best_i <- i
                best_j <- j
            }
        }
    }
    
    part = matrix(0, npoints, 1)

    switch(best_type,
        'single' = {
            part <- matrix(1, npoints, 1)
        },
        'double' = {
            part[1:best_i] <- 1
            part[(best_i + 1):length(part)] <- 2
        },
        'triple' = {
            part[1:best_i] <- 1
            part[(best_i + 1):best_j] <- 2
            part[(best_j + 1):length(part)] <- 3
        }
    )
    return(part)
}