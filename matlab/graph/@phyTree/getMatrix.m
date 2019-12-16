function [cm lab dist] = getmatrix(tr,varargin)
%GETMATRIX converts a Phytree Object into a relationship matrix.
%
%   [MATRIX, ID, DISTANCES] = GETMATRIX(T) converts the phylogenetic tree
%   object T into a logical sparse matrix, where 1's indicate that a branch
%   node (row index) is connected to its child (column index). The child
%   can be either another branch node or a leaf node. ID is a list of the
%   labels that correspond to the rows and columns of MATRIX, first the
%   leaf nodes from 1 to NUMLEAVES, then the branch nodes from NUMLEAVES+1
%   to NUMLEAVES+NUMBRANCHES, being the root the last node. DISTANCES is  
%   a column vector with one entry for every nonzero entry in MATRIX
%   traversed columnwise and representing the distance between the branch
%   node and the child. 
%
%   Example:
%
%      T = phytreeread('pf00002.tree')
%      [MATRIX ID DIST] = getmatrix(T); 
%
%  See also PHYTREE, PHYTREE/GET, PHYTREE/PDIST, PHYTREE/PRUNE, PHYTREETOOL.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/06/16 20:06:43 $


 %%% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error('Bioinfo:phytree:getmatrix:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'no_input_arguments'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
             error('Bioinfo:phytree:getmatrix:StructParamError',...
                'parameter cannot be a struct');
        end
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error('Bioinfo:phytree:getmatrix:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:phytree:getmatrix:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
            end
        end
    end
end

numLeaves = size(tr.tree,1)+1;
numNodes = numLeaves + numLeaves -1;
cm = sparse(repmat(numLeaves+1:numNodes,1,2),tr.tree(:),true,numNodes,numNodes);
if nargout>1
    lab = tr.names;
end
if nargout>2
    dist = tr.dist;
end
