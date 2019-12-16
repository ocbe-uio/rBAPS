function [ptrs,dist,names] = getcanonical(tr)
%GETCANONICAL Calculates the canonical form of a phylogenetic tree.
%
%  PTRS = GETCANONICAL(TREE) Returns the pointers of the canonical form of
%  a phylogenetic tree. In a canonical tree the leaves are ordered
%  alphabetically and the branches are ordered first by their width and
%  then alphabetically by their first element. A canonical tree is
%  isomorphic to all the trees with the same skeleton independently of the
%  order of their leaves and branches. 
%
%  [PTRS,DIST,NAMES] = GETCANONICAL(TREE) Returns also the re-ordered
%  distances and node names.
%
%  Example:
%    % create two trees with same skeleton but slightly different distances
%    b = [1 2; 3 4; 5 6; 7 8;9 10];
%    tr_1 = phytree(b,[.1 .2 .3 .3 .4 ]');
%    tr_2 = phytree(b,[.2 .1 .2 .3 .4 ]');
%    plot(tr_1)
%    plot(tr_2)
%     
%    % compare if the two trees are isomorphic
%    isequal(getcanonical(tr_1),getcanonical(tr_2))
%
%   See also PHYTREE, PHYTREEREAD, PHYTREE/GETBYNAME, PHYTREE/SELECT,
%   PHYTREE/SUBTREE.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.8.2 $ $Author: batserve $ $Date: 2006/06/16 20:06:42 $

if numel(tr)~=1
     error('Bioinfo:phytree:getcanonical:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

[dummy,h] = sort(tr.names(1:numLeaves)); %#ok
h(h)=1:numLeaves;

% compute the branch width and the first element for each one
branchWidth = ones(numLabels,1);
firstElement = [h;inf(numBranches,1)];
for ind = 1:numBranches
    branchWidth(numLeaves+ind) = sum(branchWidth(tr.tree(ind,:)));
    firstElement(numLeaves+ind) = min(firstElement(tr.tree(ind,:)));
end 

% find out how to re-order
[dummy,ord]=sortrows([branchWidth firstElement]); %#ok
iord(ord) = 1:numLabels; 

% re-order pointers
ptrs = sort(iord(tr.tree(ord(numLeaves+1:numLabels)-numLeaves,:)),2);

if nargout > 1
    dist = tr.dist(ord);
end
if nargout > 2
    names = tr.names(ord);
end


