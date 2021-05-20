clusternum <- function(X, T, k, c) {
	m <- size(X, 1) + 1
	while (!isempty(k)) {
		# Get the children of nodes at this level
		children <- X[k, 1:2]

		# Assign this node number to leaf children
		t <- (children <= m)
		T[children[t]] <- c
		# Move to next level
		k <- children[!t] - m
	}
	return(T)
}