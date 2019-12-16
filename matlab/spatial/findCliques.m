function [cliques, separators, G] = findCliques(M)
%Muuttaa graafin M kolmioituvaksi ja laskee siitä klikit ja
%separaattorit.
%Hyödynnetään Kevin Murphyn algoritmeja Graph Theory toolboxista.
%Päivitetty 12.8.2005

order=elim_order(M,ones(length(M)));
[G,cliques]=triangulate(M,order);
[jtree,root]=cliques_to_jtree(cliques,ones(length(M)));
ncliq=length(cliques);
separators=cell(ncliq-1,1);     %n-solmuisessa puussa n-1 viivaa

jono=zeros(length(ncliq));
jono(1)=root;
i=1;    
pointer=2;                      %Seuraava tyhjä paikka

while ~isempty(find(jono~=0))   %Puun leveyssuuntainen läpikäynti
    lapset=find(jtree(jono(i),:)~=0);
    jtree(:,jono(i))=0;          %Klikki käsitelty
    jono(pointer:pointer+length(lapset)-1)=lapset;
    for j=1:length(lapset)
        ehdokas = myintersect(cliques{jono(i)},cliques{lapset(j)});
        kelpaa = 1;
        for k = 1:(pointer+j-3)
            % Tutkitaan, että separaattoriehdokasta ei vielä käsitelty
            if isequal(ehdokas,separators{k})
                kelpaa = 0;
            end
        end
        if kelpaa
            separators{pointer+j-2} = ehdokas;
        end
    end
    jono(i)=0;
    pointer=pointer+length(lapset);
    i=i+1;
end

notEmpty=zeros(ncliq-1,1);
for i=1:ncliq-1
    if ~isempty(separators{i})
        notEmpty(i)=1;
    end
end
notEmpty=find(notEmpty==1);
separators=separators(notEmpty);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function order = elim_order(G, node_sizes)
% BEST_FIRST_ELIM_ORDER Greedily search for an optimal elimination order.
% order = best_first_elim_order(moral_graph, node_sizes)
%
% Find an order in which to eliminate nodes from the graph in such a way as to try and minimize the
% weight of the resulting triangulated graph.  The weight of a graph is the sum of the weights of each
% of its cliques; the weight of a clique is the product of the weights of each of its members; the
% weight of a node is the number of values it can take on.
%
% Since this is an NP-hard problem, we use the following greedy heuristic:
% at each step, eliminate that node which will result in the addition of the least
% number of fill-in edges, breaking ties by choosing the node that induces the lighest clique.
% For details, see
% - Kjaerulff, "Triangulation of graphs -- algorithms giving small total state space",
%      Univ. Aalborg tech report, 1990 (www.cs.auc.dk/~uk)
% - C. Huang and A. Darwiche, "Inference in Belief Networks: A procedural guide",
%      Intl. J. Approx. Reasoning, 11, 1994
%

% Warning: This code is pretty old and could probably be made faster.

n = length(G);
%if nargin < 3, stage = { 1:n }; end % no constraints

% For long DBNs, it may be useful to eliminate all the nodes in slice t before slice t+1.
% This will ensure that the jtree has a repeating structure (at least away from both edges).
% This is why we have stages.
% See the discussion of splicing jtrees on p68 of
% Geoff Zweig's PhD thesis, Dept. Comp. Sci., UC Berkeley, 1998.
% This constraint can increase the clique size significantly.

MG = G; % copy the original graph
uneliminated = ones(1,n);
order = zeros(1,n);
%t = 1;  % Counts which time slice we are on        
for i=1:n
  U = find(uneliminated);
  %valid = myintersect(U, stage{t});
  valid = U;
  % Choose the best node from the set of valid candidates
  min_fill = zeros(1,length(valid));
  min_weight = zeros(1,length(valid));
  for j=1:length(valid)
    k = valid(j);
    nbrs = myintersect(neighbors(G, k), U);
    l = length(nbrs);
    M = MG(nbrs,nbrs);
    min_fill(j) = l^2 - sum(M(:)); % num. added edges
    min_weight(j) = prod(node_sizes([k nbrs])); % weight of clique
  end
  lightest_nbrs = find(min_weight==min(min_weight));
  % break ties using min-fill heuristic
  best_nbr_ndx = argmin(min_fill(lightest_nbrs));
  j = lightest_nbrs(best_nbr_ndx); % we will eliminate the j'th element of valid
  %j1s = find(score1==min(score1));
  %j = j1s(argmin(score2(j1s)));
  k = valid(j);
  uneliminated(k) = 0;
  order(i) = k;
  ns = myintersect(neighbors(G, k), U);
  if ~isempty(ns)
    G(ns,ns) = 1;
    G = setdiag(G,0);
  end
  %if ~any(logical(uneliminated(stage{t}))) % are we allowed to the next slice?
  %  t = t + 1;
  %end   
end

%--------------------------------------------------------------------------

function [G, cliques, fill_ins] = triangulate(G, order)
% TRIANGULATE Ensure G is triangulated (chordal), i.e., every cycle of length > 3 has a chord.
% [G, cliques, fill_ins, cliques_containing_node] = triangulate(G, order)
% 
% cliques{i} is the i'th maximal complete subgraph of the triangulated graph.
% fill_ins(i,j) = 1 iff we add a fill-in arc between i and j.
%
% To find the maximal cliques, we save each induced cluster (created by adding connecting
% neighbors) that is not a subset of any previously saved cluster. (A cluster is a complete,
% but not necessarily maximal, set of nodes.)

MG = G;
n = length(G);
eliminated = zeros(1,n);
cliques = {};
for i=1:n
  u = order(i);
  U = find(~eliminated); % uneliminated
  nodes = myintersect(neighbors(G,u), U); % look up neighbors in the partially filled-in graph
  nodes = myunion(nodes, u); % the clique will always contain at least u
  G(nodes,nodes) = 1; % make them all connected to each other
  G = setdiag(G,0);  
  eliminated(u) = 1;
  
  exclude = 0;
  for c=1:length(cliques)
    if mysubset(nodes,cliques{c}) % not maximal
      exclude = 1;
      break;
    end
  end
  if ~exclude
    cnum = length(cliques)+1;
    cliques{cnum} = nodes;
  end
end

%fill_ins = sparse(triu(max(0, G - MG), 1));
fill_ins=1;

%--------------------------------------------------------------------------

function [jtree, root, B, w] = cliques_to_jtree(cliques, ns)
% MK_JTREE Make an optimal junction tree.
% [jtree, root, B, w] = mk_jtree(cliques, ns)
% 
% A junction tree is a tree that satisfies the jtree property, which says:
% for each pair of cliques U,V with intersection S, all cliques on the path between U and V
% contain S. (This ensures that local propagation leads to global consistency.)
%
% We can create a junction tree by computing the maximal spanning tree of the junction graph.
% (The junction graph connects all cliques, and the weight of an edge (i,j) is
% |C(i) intersect C(j)|, where C(i) is the i'th clique.)
%
% The best jtree is the maximal spanning tree which minimizes the sum of the costs on each edge,
% where cost(i,j) = w(C(i)) + w(C(j)), and w(C) is the weight of clique C,
% which is the total number of values C can take on.
%
% For details, see
%  - Jensen and Jensen, "Optimal Junction Trees", UAI 94.
%
% Input:
%  cliques{i} = nodes in clique i
%  ns(i) = number of values node i can take on
% Output:
%  jtree(i,j) = 1 iff cliques i and j aer connected
%  root = the clique that should be used as root
%  B(i,j) = 1 iff node j occurs in clique i
%  w(i) = weight of clique i



num_cliques = length(cliques);
w = zeros(num_cliques, 1); 
B = sparse(num_cliques, 1);
for i=1:num_cliques
  B(i, cliques{i}) = 1;
  w(i) = prod(ns(cliques{i}));
end


% C1(i,j) = length(intersect(cliques{i}, cliques{j})); 
% The length of the intersection of two sets is the dot product of their bit vector representation.
C1 = B*B';
C1 = setdiag(C1, 0);

% C2(i,j) = w(i) + w(j)
num_cliques = length(w);
W = repmat(w, 1, num_cliques);
C2 = W + W';
C2 = setdiag(C2, 0);

jtree = sparse(minimum_spanning_tree(-C1, C2)); % Using -C1 gives *maximum* spanning tree

% The root is arbitrary, but since the first pass is towards the root,
% we would like this to correspond to going forward in time in a DBN.
root = num_cliques;

%--------------------------------------------------------------------------


function C = myintersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = myintersect(A,B)

A = A(:)'; B = B(:)';

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 | mb==0
  C = [];
else
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  C = B(logical(bits(B)));  
end

%sum( bitget( bitand( cliquesb(i), cliquesb(j) ), 1:52 ) );

%--------------------------------------------------------------------------

function ns = neighbors(adj_mat, i)
% NEIGHBORS Find the parents and children of a node in a graph.
% ns = neighbors(adj_mat, i)

%ns = myunion(children(adj_mat, i), parents(adj_mat, i));
ns = find(adj_mat(i,:));

%--------------------------------------------------------------------------

function C = myunion(A,B)
% MYUNION Union of two sets of positive integers (much faster than built-in union)
% C = myunion(A,B)

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 & mb==0
  C = [];
elseif ma==0 & mb>0
  C = B;
elseif ma>0 & mb==0
  C = A;
else
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  bits(B) = 1;
  C = find(bits);
end

%--------------------------------------------------------------------------


function ps = parents(adj_mat, i)
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = find(adj_mat(:,i))';

%--------------------------------------------------------------------------

function cs = children(adj_mat, i, t)
% CHILDREN Return the indices of a node's children in sorted order
% c = children(adj_mat, i, t)
%
% t is an optional argument: if present, dag is assumed to be a 2-slice DBN

if nargin < 3 
  cs = find(adj_mat(i,:));
else
  if t==1
    cs = find(adj_mat(i,:));
  else
    ss = length(adj_mat)/2;
    j = i+ss;
    cs = find(adj_mat(j,:)) + (t-2)*ss;
  end
end

%--------------------------------------------------------------------------

function p=mysubset(small,large)
% MYSUBSET Is the small set of +ve integers a subset of the large set?
% p = mysubset(small, large)

% Surprisingly, this is not built-in.

if isempty(small)
  p = 1; % isempty(large);
else
  p = length(myintersect(small,large)) == length(small);
end

%--------------------------------------------------------------------------

function A = minimum_spanning_tree(C1, C2)
%
% Find the minimum spanning tree using Prim's algorithm.
% C1(i,j) is the primary cost of connecting i to j.
% C2(i,j) is the (optional) secondary cost of connecting i to j, used to break ties.
% We assume that absent edges have 0 cost.
% To find the maximum spanning tree, used -1*C.
% See Aho, Hopcroft & Ullman 1983, "Data structures and algorithms", p 237.

% Prim's is O(V^2). Kruskal's algorithm is O(E log E) and hence is more efficient
% for sparse graphs, but is implemented in terms of a priority queue.

% We partition the nodes into those in U and those not in U.
% closest(i) is the vertex in U that is closest to i in V-U.
% lowcost(i) is the cost of the edge (i, closest(i)), or infinity is i has been used.
% In Aho, they say C(i,j) should be "some appropriate large value" if the edge is missing.
% We set it to infinity.
% However, since lowcost is initialized from C, we must distinguish absent edges from used nodes.

n = length(C1);
if nargin==1, C2 = zeros(n); end
A = zeros(n);

closest = ones(1,n);
used = zeros(1,n); % contains the members of U
used(1) = 1; % start with node 1
C1(find(C1==0))=inf;
C2(find(C2==0))=inf;
lowcost1 = C1(1,:);
lowcost2 = C2(1,:);

for i=2:n
  ks = find(lowcost1==min(lowcost1));
  k = ks(argmin(lowcost2(ks)));
  A(k, closest(k)) = 1;
  A(closest(k), k) = 1;
  lowcost1(k) = inf;
  lowcost2(k) = inf;
  used(k) = 1;
  NU = find(used==0);
  for ji=1:length(NU)
    for j=NU(ji)
      if C1(k,j) < lowcost1(j)
	lowcost1(j) = C1(k,j);
	lowcost2(j) = C2(k,j);
	closest(j) = k;
      end
    end
  end
end

%--------------------------------------------------------------------------

function indices = argmin(v)
% ARGMIN Return as a subscript vector the location of the smallest element of a multidimensional array v.
% indices = argmin(v)
%
% Returns the first minimum in the case of ties.
% Example:
% X = [2 8 4; 7 3 9];
% argmin(X) = [1 1], i.e., row 1 column 1

[m i] = min(v(:));
indices = ind2subv(mysize(v), i);

%--------------------------------------------------------------------------

function M = setdiag(M, v)
% SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
% M = set_diag(M, v)

n = length(M);
if length(v)==1
  v = repmat(v, 1, n);
end

% e.g., for 3x3 matrix,  elements are numbered
% 1 4 7 
% 2 5 8 
% 3 6 9
% so diagnoal = [1 5 9]


J = 1:n+1:n^2;
M(J) = v;

%-------------------------------------------------------------------------

function sz = mysize(M)
% MYSIZE Like the built-in size, except it returns n if M is a vector of length n, and 1 if M is a scalar.
% sz = mysize(M)
% 
% The behavior is best explained by examples
% - M = rand(1,1),   mysize(M) = 1,      size(M) = [1 1]
% - M = rand(2,1),   mysize(M) = 2,      size(M) = [2 1]
% - M = rand(1,2),   mysize(M) = 2,      size(M) = [1 2]
% - M = rand(2,2,1), mysize(M) = [2 2],  size(M) = [2 2]
% - M = rand(1,2,1), mysize(M) = 2,      size(M) = [1 2]

if myisvector(M)
  sz = length(M);
else
  sz = size(M);
end

%--------------------------------------------------------------------------

function sub = ind2subv(siz, ndx)
% IND2SUBV Like the built-in ind2sub, but returns the answer as a row vector.
% sub = ind2subv(siz, ndx)
%
% siz and ndx can be row or column vectors.
% sub will be of size length(ndx) * length(siz).
%
% Example
% ind2subv([2 2 2], 1:8) returns
%  [1 1 1
%   2 1 1
%   ...
%   2 2 2]
% That is, the leftmost digit toggle fastest.
%
% See also SUBV2IND

n = length(siz);

if n==0
  sub = ndx;
  return;
end  

if all(siz==2)
  sub = dec2bitv(ndx-1, n);
  sub = sub(:,n:-1:1)+1;
  return;
end

cp = [1 cumprod(siz(:)')];
ndx = ndx(:) - 1;
sub = zeros(length(ndx), n);
for i = n:-1:1 % i'th digit
  sub(:,i) = floor(ndx/cp(i))+1;
  ndx = rem(ndx,cp(i));
end

%%%%%%%%%%

function bits = dec2bitv(d,n)
% DEC2BITV Convert a decimal integer to a bit vector.
% bits = dec2bitv(d,n) is just like the built-in dec2bin, except the answer is a vector, not a string.
% n is an optional minimum length on the bit vector.
% If d is a vector,  each row of the output array will be a bit vector.


if (nargin<2)
  n=1; % Need at least one digit even for 0.
end
d = d(:);

[f,e]=log2(max(d)); % How many digits do we need to represent the numbers?
bits=rem(floor(d*pow2(1-max(n,e):0)),2);


%------------------------------------------------------------------------

function r = myisvector(V)
%Kuten isvector(V)

A = size(V);
r = (length(A) == 2) & (min(A) == 1);
