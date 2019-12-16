function tr = prettyOrder(tr)
%PRETTYORDER Reorders the leaf nodes to avoid branch crossings.
%
%   T2 = PRETTYORDER(T1) Reorders the leaf nodes in the phylogenetic tree
%   T1 such that the layout of the tree does not contain branch crossings.

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.8.1 $ $Author: batserve $ $Date: 2005/06/09 21:56:11 $

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

L = [ones(numLeaves,1); zeros(numBranches,1)];
for ind = 1 : numBranches
   L(ind+numLeaves) = sum(L(tr.tree(ind,:)));
end
X = zeros(numLabels,1);
for ind = numBranches:-1:1
   X(tr.tree(ind,:)) = tr.dist(tr.tree(ind,:))+X(ind+numLeaves);
end
Li = zeros(1,numLabels); Ls = Li;
Ls(numLabels) = numLeaves; 
for ind = numBranches:-1:1
    Ls(tr.tree(ind,:)) = Ls(ind+numLeaves);
    Li(tr.tree(ind,:)) = Li(ind+numLeaves);
    if diff(X(tr.tree(ind,:)))>=0
        Ls(tr.tree(ind,1)) = Li(tr.tree(ind,1)) + L(tr.tree(ind,1));
        Li(tr.tree(ind,2)) = Ls(tr.tree(ind,2)) - L(tr.tree(ind,2));
    else
        Ls(tr.tree(ind,2)) = Li(tr.tree(ind,2)) + L(tr.tree(ind,2));
        Li(tr.tree(ind,1)) = Ls(tr.tree(ind,1)) - L(tr.tree(ind,1));
    end
end
  
tr.names(Ls(1:numLeaves))=tr.names(1:numLeaves);
tr.dist(Ls(1:numLeaves))=tr.dist(1:numLeaves);
Ls(numLeaves+1:numLabels)=numLeaves+1:numLabels;
tr.tree = Ls(tr.tree);
