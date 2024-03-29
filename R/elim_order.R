elim_order <- function(G, node_sizes) {
  #  BEST_FIRST_ELIM_ORDER Greedily search for an optimal elimination order.
  #  order <- best_first_elim_order(moral_graph, node_sizes)
  #  Find an order in which to eliminate nodes from the graph in such a way as to try and minimize the
  #  weight of the resulting triangulated graph.  The weight of a graph is the sum of the weights of each
  #  of its cliques the weight of a clique is the product of the weights of each of its members the
  #  weight of a node is the number of values it can take on.
  #  Since this is an NP - hard problem, we use the following greedy heuristic:
  #  at each step, eliminate that node which will result in the addition of the least
  #  number of fill - in edges, breaking ties by choosing the node that induces the lighest clique.
  #  For details, see
  #  - Kjaerulff, "Triangulation of graphs  - - algorithms giving small total state space",
  #       Univ. Aalborg tech report, 1990 (www.cs.auc.dk/!uk)
  #  - C. Huang and A. Darwiche, "Inference in Belief Networks: A procedural guide",
  #       Intl. J. Approx. Reasoning, 11, 1994

  #  Warning: This code is pretty old and could probably be made faster.

  n <- length(G)
  # if (nargin < 3, stage = { 1:n } end# no constraints) {

  #  For long DBNs, it may be useful to eliminate all the nodes in slice t before slice t + 1.
  #  This will ensure that the jtree has a repeating structure (at least away from both edges).
  #  This is why we have stages.
  #  See the discussion of splicing jtrees on p68 of
  #  Geoff Zweig's PhD thesis, Dept. Comp. Sci., UC Berkeley, 1998.
  #  This constraint can increase the clique size significantly.

  MG <- G# copy the original graph
  uneliminated <- ones(1, n)
  order <- zeros(1, n)
  # t <- 1 # Counts which time slice we are on
  for (i in 1:n) {
    U <- find(uneliminated)
    # valid <- myintersect(U, stage{t})
    valid <- U
    # Choose the best node from the set of valid candidates
    min_fill <- zeros(1, length(valid))
    min_weight <- zeros(1, length(valid))
    for (j in 1:length(valid)) {
      k <- valid(j)
      nbrs <- myintersect(neighbors(G, k), U)
      l <- length(nbrs)
      M <- MG[nbrs, nbrs]
      min_fill[j] <- l^2 - sum(M)# num. added edges
      min_weight[j] <- prod(node_sizes[k, nbrs])# weight of clique
    }
    lightest_nbrs <- find(min_weight == min(min_weight))
    # break ties using min - fill heuristic
    best_nbr_ndx <- argmin(min_fill[lightest_nbrs])
    j <- lightest_nbrs[best_nbr_ndx] # we will eliminate the j'th element of valid
    # j1s <- find(score1 == min(score1))
    # j <- j1s(argmin(score2(j1s)))
    k <- valid(j)
    uneliminated[k] <- 0
    order[i] <- k
    ns <- myintersect(neighbors(G, k), U)
    if (!is.null(ns)) {
      G[ns, ns] <- 1
      G <- setdiag(G, 0)
    }
    # if (!any(as.logical(uneliminated(stage{t})))# are we allowed to the next slice?) {
    #   t <- t + 1
    # }
  }
  return(order)
}
