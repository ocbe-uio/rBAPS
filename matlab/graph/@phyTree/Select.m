function [sel,sell,selb] = select(tr,varargin)
%SELECT Selects tree branches and leaves.
%
%   S = SELECT(T,N) returns a logical vector S of size [NUMNODES x 1]
%   indicating N closest nodes to the root node of the phylogenetic tree
%   object T (NUMNODES = NUMLEAVES + NUMBRANCHES). The first criteria used
%   is branch levels, then patristic distance, also known as tree distance.
%   By default SELECT, uses INF as the value of N, therefore SELECT(T) will
%   return a vector of 'trues'.
%
%   S = SELECT(...,'REFERENCE',R) changes the reference point(s) to measure
%   the closeness. R can be 'root' (default) or 'leaves'. When using
%   'leaves', a node to be tested may have different distances to its
%   descendant leaves, which are the references (e.g. a non-ultrametric
%   tree), if this the case the minimum distance to any descendant leaf
%   will be considered. R may also be an index which points to any node of
%   the tree.
%
%   S = SELECT(...,'CRITERIA',C) changes the criteria used to measure
%   closeness. If C='levels' (default) then the first criteria is branch
%   levels and then patristic distance. If C='distance' then the first
%   criteria is patristic distance and then branch levels.
%
%   S = SELECT(...,'THRESHOLD',V) selects all the nodes which closeness is
%   less or equal than the threshold value V. Observe that either
%   'criteria' or either 'reference' can be used. If N is not specified
%   N = INF, otherwise the output can be further size limited by N.
%
%   S = SELECT(...,'EXCLUDE',E) sets a post-filter which excludes all the
%   branch nodes from S when E=='branches' or all the leaf nodes when
%   E=='leaves'. The default is 'none'.
%
%   S = SELECT(...,'PROPAGATE',P) activates a post-functionality which
%   propagates the selected nodes to the leaves when P=='toleaves' or
%   towards the root finding a common ancestor when P=='toroot'. The
%   default is 'none', P may also be 'both'. 'PROPAGATE' switch acts after
%   'EXCLUDE' switch.  
%
%   [S,SELLEAVES,SELBRANCHES] = SELECT(...) returns two additional logical
%   vectors, one for the selected leaves and one for the selected branches.
%
%  Examples:
%
%      % Load a phylogenetic tree created from a protein family:
%      tr = phytreeread('pf00002.tree');
%  
%      % To find close products for a given protein (e.g. vips_human):
%      ind = getbyname(tr,'vips_human');
%      [sel,sel_leaves] = select(tr,'criteria','distance',...
%                                'threshold',0.6,'reference',ind);
%      view(tr,sel_leaves)
%
%      % To find potential outliers in the tree use:
%      [sel,sel_leaves] = select(tr,'criteria','distance','threshold',.3,...
%           'reference','leaves','exclude','leaves','propagate','toleaves');
%      view(tr,~sel_leaves)
%
%  
%  See also PHYTREE, PHYTREE/GET, PHYTREE/PDIST, PHYTREE/PRUNE, PHYTREETOOL.

%
% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.6.7.2.1 $ $Author: batserve $ $Date: 2006/07/27 21:37:50 $

if numel(tr)~=1
     error('Bioinfo:phytree:select:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

% set defaults
V=inf;
CriteriaIsDistance = false;
ReferenceIs = 'root';
ExcludeSwitch = false;
PostPropagate = false;

% check is first argument is N, otherwise set N to default
if (nargin>1 && isnumeric(varargin{1}) && isreal(varargin{1}))
    N = floor(varargin{1});
    first_arg = 3;
else
    N = inf;
    first_arg = 2;
end

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% identify input arguments
if nargin - first_arg + 1 > 0
   if rem(nargin - first_arg,2) == 0
       error('Bioinfo:phytree:select:IncorrectNumberOfArguments',...
             'Incorrect number of arguments to %s.',mfilename);
   end
   okargs = {'reference','criteria','threshold','exclude','propagate'};
   for j=first_arg - 1 : 2 : nargin - first_arg + 1
       pname = varargin{j};
       pval = varargin{j+1};
       k = find(strncmpi(pname,okargs,numel(pname)));
       if isempty(k)
           error('Bioinfo:phytree:select:UnknownParameterName',...
               'Unknown parameter name: %s.',pname);
       elseif length(k)>1
           error('Bioinfo:phytree:select:AmbiguousParameterName',...
               'Ambiguous parameter name: %s.',pname);
       else
           switch(k)
               case 1  % reference
                   if islogical(pval)
                       if numel(pval)==numLabels 
                          pval = pval(:)==true;
                       elseif numel(pval)==numLeaves
                          pval = [pval(:);false(numBranches,1)];
                       elseif numel(pval)==numBranches
                          pval = [false(numLeaves,1);pval(:)];
                       else
                         error('Bioinfo:phytree:select:InvalidSizeLogicalReferenceNode',...
                               'When reference node is a logical vector it must contain NUMNODES, NUMLEAVES or NUMBRANCHES elements.')
                       end
                       pval = find(pval);
                       ReferenceIs = 'node';
                       if numel(pval) ~= 1 
                           error('Bioinfo:phytree:select:InvalidValuesLogicalReferenceNode',...
                               'When reference node is a logical vector one element must be true and all others false.')
                       else
                           ReferenceNode = pval;
                       end                       
                   elseif isnumeric(pval)
                       ReferenceIs = 'node';
                       if numel(pval) ~= 1 
                           error('Bioinfo:phytree:select:InvalidSizeReferenceNode',...
                               'Reference node must be scalar.')
                       elseif all(pval~=1:numLabels)
                           error('Bioinfo:phytree:select:InvalidValueReferenceNode',...
                               'Incorrect reference node.')
                       else
                           ReferenceNode = pval;
                       end
                   else
                       h = strmatch(lower(pval),{'root','leaves'}); %#ok
                       if numel(h) 
                           switch(h)
                               case 1
                                   ReferenceIs = 'root';
                               case 2
                                   ReferenceIs = 'leaves';
                           end
                       else error('Bioinfo:phytree:select:InvalidStringReferenceNode',...
                             'Invalid string for the reference node.');
                       end
                   end
                case 2  % criteria
                   h = strmatch(lower(pval),{'distance','levels'}); %#ok
                   if numel(h) 
                       CriteriaIsDistance = (h == 1);
                   else error('Bioinfo:phytree:select:InvalidCriteria',...
                              'Invalid string for criteria.');
                   end
               case 3 % threshold
                   V = pval;
                   if (~isnumeric(V) || ~isreal(V) || numel(V)>1)
                       error('Bioinfo:phytree:select:InvalidThreshold',...
                             'Invalid value for V.');  
                   end
               case 4 % exclude
                   h = strmatch(lower(pval),{'branches','leaves','none'}); %#ok
                   if numel(h)
                       switch(h)
                           case 1
                               ExcludeSwitch = true;
                               ExcludeType = 'branches';
                           case 2
                               ExcludeSwitch = true;
                               ExcludeType = 'leaves';
                           case 3
                               ExcludeSwitch = false;
                       end
                       else error('Bioinfo:phytree:select:InvalidExcludeOption',...
                              'Invalid string for exclude switch.');
                   end
               case 5 % propagate
                   h = strmatch(lower(pval),{'toleaves','toroot','both','none'}); %#ok
                   if numel(h)
                       switch(h)
                           case 1
                               PostPropagate = true;
                               PostPropagateType = 'toleaves';
                           case 2
                               PostPropagate = true;
                               PostPropagateType = 'toroot';
                           case 3
                               PostPropagate = true;
                               PostPropagateType = 'both';
                           case 4
                               PostPropagate = false;
                       end
                       else error('Bioinfo:phytree:select:InvalidPostPropagate',...
                              'Invalid string for post-propagate switch.');
                   end
           end % switch(k)
       end %  if ...
   end % for j=...
end % nargin  

% calculate the distance (and levels) of every node to the reference
levels2Ref = zeros(numLabels,1);
switch ReferenceIs
    case 'root'
        % calculate the distance to the root for every node
        dist2Ref = tr.dist;
        for ind = numBranches:-1:1
            dist2Ref(tr.tree(ind,:)) = ...
                dist2Ref(tr.tree(ind,:)) + dist2Ref(ind+numLeaves);
            levels2Ref(tr.tree(ind,:)) = levels2Ref(ind+numLeaves) + 1;
        end
    case 'leaves'
        % calculate the distance to the closest leaf for every node
        dist2Ref = zeros(numLabels,1);
        for ind = 1:numBranches
            dist2Ref(ind+numLeaves) = ...
                min(dist2Ref(tr.tree(ind,:))+tr.dist(tr.tree(ind,:)));
            levels2Ref(ind+numLeaves) = min(levels2Ref(tr.tree(ind,:))) + 1;
        end
    case 'node'
        refVector = zeros(numLabels,1);
        refVector(ReferenceNode)=1;
        dist2Ref = pdist(tr,'SquareForm',true,'nodes','all')...
                   * refVector;
        tr.dist = ones(numLabels,1); % to count now levels !
        levels2Ref = pdist(tr,'SquareForm',true,'nodes','all') ...
                   * refVector;
end

% applies the threshold value
if CriteriaIsDistance
    sel = dist2Ref < V;
else % ~CriteriaIsDistance
    sel = levels2Ref < V;
end % if CriteriaIsDistance

% needs to remove additional nodes because of N
if sum(sel)>N
    if CriteriaIsDistance
        [dum,h]=sortrows([dist2Ref levels2Ref]); %#ok
    else % ~CriteriaIsDistance
        [dum,h]=sortrows([levels2Ref dist2Ref]); %#ok
    end % if CriteriaIsDistance
    g=h(sel(h));
    sel(g(N+1:end))=false;
end

% exclude option
if ExcludeSwitch
    switch ExcludeType
        case 'branches'
            sel((1+numLeaves):numLabels) = false;
        case 'leaves'
            sel(1:numLeaves) = false;
    end
end

% post-propagate option
if PostPropagate

    % expands all the current nodes towards the leaves
    if any(strcmp({'toleaves','both'},PostPropagateType))
        for ind = numBranches:-1:1
            if sel(ind+numLeaves)
                sel(tr.tree(ind,:))=true;
            end
        end
    end

    % propagates towards the root finding the common ancestors
    if any(strcmp({'toroot','both'},PostPropagateType))
        
        % find closest common branch for every pair of nodes
        % diagonal is invalid ! but not needed

        % initializing full matrix
        commf = zeros(numLabels,'int16');
        children = false(1,numLabels);
        for ind = numBranches:-1:1
            children(:) = false;
            children(ind+numLeaves) = true;
            for ind2 = ind:-1:1
                if children(ind2+numLeaves)
                    children(tr.tree(ind2,:))=true;
                end
            end
            commf(children,children)=int16(ind);
        end
        commf = commf(sel,sel);
        commf = commf - diag(diag(commf));
        commf = unique(commf(commf(:)>0));
        sel(commf+numLeaves) = true;
        
        % now propagates towards the common ancestor
        for ind = 1:double(max(commf))
            if any(sel(tr.tree(ind,:)))
                sel(ind+numLeaves) = true;
            end
        end
    end
end

if nargout > 1
    sell = sel(1:numLeaves);
end
if nargout > 2
    selb = sel(numLeaves+1:numLabels);
end
