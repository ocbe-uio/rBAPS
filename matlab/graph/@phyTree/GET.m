function varargout = get(tr,varargin)
%GET Get information about a phylogenetic tree object.
%   [VALUE1,VALUE2, ...] = GET(TREE,'NAME1','NAME2', ...) returns the
%   contents of the specified fields for the PHYTREE object TREE.
%
%   The valid choices for 'NAME' are:
%     'POINTERS'    : branch to leaf/branch connectivity list
%     'DISTANCES'   : edge length for every leaf/branch
%     'NUMLEAVES'   : number of leaves
%     'NUMBRANCHES' : number of branches
%     'NUMNODES'    : number of nodes (numleaves + numbranches)
%     'LEAFNAMES'   : names of the leaves
%     'BRANCHNAMES' : names of the branches
%     'NODENAMES'   : names of all the nodes
%
%   GET(TREE) displays all property names and their current values for
%   the PHYTREE object TREE.
% 
%   V = GET(TREE) returns a structure where each field name is the name of
%   a property of TREE and each field contains the value of that property.
%
%   Examples:
%     tr = phytreeread('pf00002.tree')
%     protein_names = get(tr,'LeafNames')
%
%   See also PHYTREE, PHYTREEREAD, PHYTREE/SELECT, PHYTREE/GETBYNAME.

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.6.7 $ $Author: batserve $ $Date: 2005/06/09 21:55:54 $

if numel(tr)~=1
     error('Bioinfo:phytree:get:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves;

% get input without arguments displays a summary
if nargin == 1
   s.NumLeaves = numLeaves;
   s.NumBranches = numBranches;
   s.NumNodes = numLabels;
   s.Pointers = tr.tree;
   s.Distances = tr.dist;
   s.LeafNames = tr.names(1:numLeaves);
   s.BranchNames = tr.names(numLeaves+1:numLabels);
   s.NodeNames = tr.names;
   if nargout == 0
       disp(s)
   else
       varargout{1} = s;
   end
   return;
end

okargs = {'pointers','distances','numleaves','numbranches',...
    'numnodes','leafnames','branchnames','nodenames'};
for ind = 2 : nargin
    pname = varargin{ind-1};
    k = find(strncmpi(pname,okargs,numel(pname)));
    if isempty(k)
        error('Bioinfo:phytree:get:UnknownParameterName',...
            'Unknown parameter name: %s.',pname);
    elseif length(k)>1
        error('Bioinfo:phytree:get:AmbiguousParameterName',...
            'Ambiguous parameter name: %s.',pname);
    else
        switch(k)
            case 1 % pointers
                varargout{ind-1} = tr.tree; %#ok
            case 2 % distances
                varargout{ind-1} = tr.dist;
            case 3 % numleaves
                varargout{ind-1} = numLeaves;
            case 4 % numbranches
                varargout{ind-1} = numBranches;
            case 5 % numNodes
                varargout{ind-1} = numLabels;
            case 6 % leafnames
                varargout{ind-1} = tr.names(1:numLeaves);
            case 7 % branchnames
                varargout{ind-1} = tr.names(numLeaves+1:numLabels);
            case 8 % nodenames
                varargout{ind-1} = tr.names;
        end
    end
end
