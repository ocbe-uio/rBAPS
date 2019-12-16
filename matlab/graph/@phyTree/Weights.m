function W = weights(tr)
%WEIGHTS Tree based sequence weights.
%
%  W = WEIGHTS(T) Calculates branch proportional weights for every leaf in
%  the tree using the Thompson-Higgins-Gibson method. The distance of every
%  segment of the tree is adjusted by dividing it by the number of leaves
%  it contains. The sequence weights are the result of normalizing to the
%  unity the new patristic distances between every leaf and the root.
%
%  Example:
%
%     % Create an ultrametric tree with specified branch distances
%     bd = [1 2 3]';
%     tr_1 = phytree([1 2;3 4;5 6],bd)
%     view(tr_1) 
%     weights(tr_1)
%
%  See also MULTIALIGN, PHYTREE, PROFALIGN, SEQLINKAGE.

% References:
%   J.D. Thompson, D.G. Higgins, and T.J. Gibson. Nucleic Acids Res. (1994)
%   22(22):4673-4680.
%   S.Henikoff and J. G. Henikoff. JMB. (1994) 243(4):574--578.
%
% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.8.2 $Author: batserve $ $Date: 2005/06/17 20:19:24 $

if numel(tr)~=1
     error('Bioinfo:phytree:weights:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% calculate the branch width
branchWidth = ones(numLabels,1);
for ind = 1:numBranches
    branchWidth(numLeaves+ind) = sum(branchWidth(tr.tree(ind,:)));
end 

% adjust the distances
tr.dist = tr.dist ./ branchWidth;

% calculate distance of every leave to root
cdist = tr.dist; 
for ind = numBranches:-1:1
    cdist(tr.tree(ind,:)) = cdist(tr.tree(ind,:)) + cdist(ind+numLeaves);
end 

W = cdist(1:numLeaves);
W = W./max(W);

