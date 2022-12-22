triangulate <- function(G, order) {
  #  TRIANGULATE Ensure G is triangulated (chordal), i.e., every cycle of length > 3 has a chord.
  #  [G, cliques, fill_ins, cliques_containing_node] = triangulate(G, order)
  #  cliques{i} is the i'th maximal complete subgraph of the triangulated graph.
  #  fill_ins[i, j] <- 1 iff we add a fill - in arc between i and j.
  #  To find the maximal cliques, we save each induced cluster (created by adding connecting
  #  neighbors) that is not a subset of any previously saved cluster. (A cluster is a complete,
  #  but not necessarily maximal, set of nodes.)

  MG <- G
  n <- length(G)
  eliminated <- zeros(1, n)
  cliques = list()
  for (i in 1:n) {
    u <- order[i]
    U <- find(!eliminated)# uneliminated
    nodes <- myintersect(neighbors(G, u), U)# look up neighbors in the partially filled - in graph   # TODO: translate neighbors
    nodes <- myunion(nodes, u)# the clique will always contain at least u # TODO: translate myunion
  ----------------------- line 21 -----------------------
    G(nodes,nodes) = 1; % make them all connected to each other
    G[nodes, nodes] <- 1# make them all connected to each other
  ----------------------- line 22 -----------------------
    G = setdiag(G,0);
    G <- setdiag(G, 0)
  ----------------------- line 23 -----------------------
    eliminated(u) = 1;
    eliminated[u] <- 1
  ----------------------- line 24 -----------------------


  ----------------------- line 25 -----------------------
    exclude = 0;
    exclude <- 0
  ----------------------- line 26 -----------------------
    for c=1:length(cliques)
    for (c in 1:length(cliques)) {
  ----------------------- line 27 -----------------------
    if mysubset(nodes,cliques{c}) % not maximal
    if (mysubset(nodes, cliques{c})# not maximal) {
  ----------------------- line 28 -----------------------
    exclude = 1;
    exclude <- 1
  ----------------------- line 29 -----------------------
    break;
    break
  ----------------------- line 30 -----------------------
    end
    }
  ----------------------- line 31 -----------------------
  end
  }
  ----------------------- line 32 -----------------------
    if ~exclude
    if (!exclude) {
  ----------------------- line 33 -----------------------
      cnum = length(cliques)+1;
      cnum <- length(cliques) + 1
  ----------------------- line 34 -----------------------
      cliques{cnum} = nodes;
      cliques{cnum} = nodes
  ----------------------- line 35 -----------------------
    end
    }
  ----------------------- line 36 -----------------------
  end
  }
  ----------------------- line 37 -----------------------


  ----------------------- line 38 -----------------------
  %fill_ins = sparse(triu(max(0, G - MG), 1));
  # fill_ins <- sparse(triu(max(0, G - MG), 1))
  ----------------------- line 39 -----------------------
  fill_ins=1;
  fill_ins=1
  ----------------------- line 40 -----------------------
  NA
  return(list("G" = G, "cliques" = cliques, "fill_ins" = fill_ins))
}
