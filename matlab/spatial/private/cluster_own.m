function T = cluster_own(Z,nclust)
% search down the dendogram from the root, until nclust clusters are found
% comments added by Lu Cheng
% 04.01.2011

maxclust = nclust;
% Start of algorithm
m = size(Z,1)+1;
T = zeros(m,1);

% maximum number of clusters based on inconsistency
if m <= maxclust
  T = (1:m)';
elseif maxclust==1
  T = ones(m,1);
else
  clsnum = 1;
  for k = (m-maxclust+1):(m-1)
     i = Z(k,1); % left tree
     if i <= m % original node, no leafs
        T(i) = clsnum;
        clsnum = clsnum + 1;
     elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
        T = clusternum(Z, T, i-m, clsnum);
        clsnum = clsnum + 1;
     end
     i = Z(k,2); % right tree
     if i <= m  % original node, no leafs
        T(i) = clsnum;
        clsnum = clsnum + 1;
     elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
        T = clusternum(Z, T, i-m, clsnum);
        clsnum = clsnum + 1;
     end
  end
end
   
function T = clusternum(X, T, k, c)
m = size(X,1)+1;
while(~isempty(k))
   % Get the children of nodes at this level
   children = X(k,1:2);
   children = children(:);

   % Assign this node number to leaf children
   t = (children<=m);
   T(children(t)) = c;
   
   % Move to next level
   k = children(~t) - m;
end

