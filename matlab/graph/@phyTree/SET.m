function tr = set(tr,varargin) %#ok
%SET  Set object properties of a phylogenetic tree object.
%
%  Properties in a phylogenetic tree object can not be manually set.
%  A PHYTREE object must be created by its constructor method PHYTREE
%  or by using one of the functions: PHYTREEREAD, SEQLINKAGE, SEQNEIGHJOIN.
%
%  See also: PHYTREE, PHYTREEREAD, SEQLINKAGE, SEQNEIGHJOIN.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.8.2.2.1 $ $Author: batserve $ $Date: 2006/07/27 21:37:51 $

% error('Bioinfo:phytree:set:NotAllowedMethod',...
%   ['Properties in a phylogenetic tree object can not be manually set.\n'...
%    'A PHYTREE object must be created by its constructor method PHYTREE\n'...
%    'or by using one of the functions: PHYTREEREAD, SEQLINKAGE, SEQNEIGHJOIN.'])

% numBranches = size(tr.tree,1);
% numLeaves = numBranches + 1;
% 
% tr.names(1:numLeaves)= varargin{2}';