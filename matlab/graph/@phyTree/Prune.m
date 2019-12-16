function tr = prune(tr,sel,varargin)
%PRUNE Reduces a phylogenetic tree by removing branch and leaf nodes.  
%   T2 = PRUNE(T1,NODES) prunes the nodes listed in the NODES vector from
%   the tree T1.  Any branch (or leaf) node listed in NODES and all their
%   descendants will disappear. The respective 'parent' nodes will be
%   connected to the respective 'brother' nodes as required. NODES in the
%   tree are indexed as [1:NUMLEAVES] for the leaves and as
%   [NUMLEAVES+1 : NUMLEAVES+NUMBRANCHES] for the branches. NODES can also
%   be a logical array of following sizes: [NUMLEAVES+NUMBRANCHES x 1],
%   [NUMLEAVES x 1] or [NUMBRANCHES x 1].
%
%   T2 = PRUNE(T1,NODES,'MODE','EXCLUSIVE') changes the pruning mode to
%   'EXCLUSIVE', i.e. only the descendants of NODES will be pruned. Then
%   NODES will become leaves as long as they do not have a predecessor in
%   the list NODES. In this case pruning is the process of reducing a tree
%   by turning some branch nodes into leaf nodes, and removing the leaf
%   nodes under the original branch. Default is 'INCLUSIVE' and it behaves
%   as explained above, i.e. the listed NODES are also pruned.
%  
%   Examples:
%
%      % Load a phylogenetic tree created from a protein family:
%      tr = phytreeread('pf00002.tree');
%      view(tr)
%      
%      % To remove all the 'mouse' proteins use:
%      ind   = getbyname(tr,'mouse');
%      tr    = prune(tr,ind);
%      view(tr)
%  
%      % To remove potential outliers in the tree use:
%      [sel,sel_leaves] = select(tr,'criteria','distance','threshold',.3,...
%           'reference','leaves','exclude','leaves','propagate','toleaves');
%      tr = prune(tr,~sel_leaves)
%      view(tr)
%      
%   See also PHYTREE, PHYTREE/SELECT, PHYTREE/GET, PHYTREETOOL.

%   Copyright 2003-2005 The MathWorks, Inc.
%   $Revision: 1.1.4.8.2.1 $    $Date: 2004/11/30 03:45:24 $

% set default
exclusiveMode = false;

if numel(tr)~=1
     error('Bioinfo:phytree:prune:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

btr = tr;
numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

% validate sel
if islogical(sel)
    if numel(sel)==numLabels 
        sel = sel(:)==true;
    elseif numel(sel)==numLeaves
        sel = [sel(:);false(numBranches,1)];
    elseif numel(sel)==numBranches
        sel = [false(numLeaves,1);sel(:)];
    else
        error('Bioinfo:IncorrectLogical',...
        'Logical vector must have the same number of elements as nodes in the Phylogenetic Tree');
    end
elseif isnumeric(sel) && isreal(sel) && all(sel>=1) && all(sel<=numLabels)
    tem(numLabels)=false;
    tem(floor(sel))=true;
    sel=tem(:);
else
    error('Bioinfo:IncorrectTypeofArguments','Invalid value for NODES');
end

nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2)
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'mode',''};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % classifiers
                    modeOptions = {'exclusive','inclusive'};
                    modeSelected = strmatch(lower(pval),modeOptions); %#ok
                    if isempty(modeSelected) 
                      error('Bioinfo:NotValidMode','Not a valid mode.')
                    end
                    exclusiveMode = modeSelected==1;
            end
        end
    end
end

% shortcut for an empty sel
if ~sum(sel) 
    return; 
end

% when inclusiveMode if the two chidren of a branch are selected then the
% parent node should also be selected
if ~exclusiveMode
    for ind = 1:numBranches
        if all(sel(tr.tree(ind,:)))
            sel(ind+numLeaves) = true;
        end
    end
end     

% find descendants not selected under selected nodes
for ind = numBranches:-1:1
    if sel(ind+numLeaves)
        sel(tr.tree(ind,:))=true;
    end
end

 if sel(numLabels) 
     warning('Bioinfo:PrunedRoot',...
             'Can not prune the root node in a Phylogenetic Tree.')
     tr=btr; return
 end

% obtain parents for every node
parents(tr.tree(:)) = repmat(numLeaves+1:numLabels,2,1)';

if ~exclusiveMode % (the selected nodes are deleted with their descendants)
    % find the top selected nodes in order to edit branches
    htop = find(~[sel(parents);0]&sel);
    % for every top node do the junction
    for ind = 1:length(htop)
        g=htop(ind);
        mypar = parents(g);
            if mypar < numLabels  % my parent is NOT the root
                                  % then connect brother to granparent
              mygrpar = parents(mypar);                       % grandparent
              myuncle = setxor(tr.tree(mygrpar-numLeaves,:),mypar); % uncle 
              mybro = setxor(tr.tree(mypar-numLeaves,:),g);       % brother
              tr.tree(mygrpar-numLeaves,:) = [myuncle mybro];
              tr.dist(mybro) = tr.dist(mybro) + tr.dist(mypar);
              parents(mybro) = mygrpar;
            end
        sel(mypar) = true; %also delete my par   
        
    end 
    if sum(~sel) == 1
            warning('Bioinfo:NotAMinimumTree',...
            'The selected nodes lead to only one leaf, Phylogenetic Tree not pruned.')
            tr=btr; return
    end
    % find indexes to change tree 
    permuta = 1:numLabels;
    permuta(sel) = [];
    ipermuta(permuta) = 1:length(permuta);
    permutaBranches = permuta(permuta>numLeaves)-numLeaves;
    % update all tree structure fields
    tr.names = tr.names(permuta);
    tr.dist = tr.dist(permuta);
    tr.dist(end) = 0;
    tr.tree = tr.tree(permutaBranches,:);
    tr.tree = ipermuta(tr.tree);
    
else % exclusiveMode (the selected nodes are not deleted, only their descendants)

    % unselect leaves which are already in the top
    sel(1:numLeaves)=sel(parents(1:numLeaves));
    % find the top selected nodes in order to edit branches
    top = [~sel(parents);1] & sel;  
    % find the new leaves (no deleted leaves + branches that become leaves)
    newLeaves = [~sel(1:numLeaves);top(numLeaves+1:end)];
    % find which branches will stay
    stayingBranches = ~sel(numLeaves + 1 : numLabels);
    % setting new indexes to change the tree architecture
    permuta =  [find(newLeaves);numLeaves+find(stayingBranches)];
    ipermuta(permuta) = 1:length(permuta);
    % update all tree structure fields
    tr.names = tr.names(permuta);
    tr.dist = tr.dist(permuta);
    tr.dist(end) = 0;
    tr.tree = tr.tree(stayingBranches,:);
    tr.tree = ipermuta(tr.tree);
    % calling phytree with this format to force edge-crossing check
    tr = phytree(tr.tree,tr.dist,tr.names);

end % if ~exclusiveMode
        

