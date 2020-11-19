clusternum <- function(X, T, k, c) {
	m <- size(X, 1) + 1
	while (!is.null(k)) {
		# % Get the children of nodes at this level
		children <- X[k, 1:2]
		children <- children[, ]

		# % Assign this node number to leaf children
		t <- (children <= m)
		T[children(t)] <- c

		# % Move to next level
		k <- children(!t) - m
	}
	return(T)
}
