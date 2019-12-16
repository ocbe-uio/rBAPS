function tr = reRoot(tr,node,distance)
%REROOT changes the root of a phylogenetic tree.
%
%   T2 = REROOT(T1) changes the root of the phylogenetic tree T1 using the
%   mid-point method. The mid-point is the location where the means of
%   the branch lengths of either side of the tree are equalized. The
%   original root is deleted.
%
%   T2 = REROOT(T1,NODE) changes the root to the branch indexed by NODE.
%   The new root is placed at half the distance between NODE and its
%   parent. 
%
%   T2 = REROOT(T1,NODE,DISTANCE) re-roots T1 by placing the new root at a
%   given DISTANCE from the reference NODE towards the root of the tree.
%
%   Note: The new branch in T2 representing the root is labeled as 'Root'.
%
%   Example:
%
%      % Create an ultrametric tree
%      tr_1 = phytree([5 7;8 9;6 11; 1 2;3 4;10 12;14 16;15 17;13 18])
%      plot(tr_1,'branchlabels',true) 
%
%      % Place the new root at 'Branch 7'
%      sel = getbyname(tr_1,'Branch 7');
%      tr_2 =  reroot(tr_1,sel)
%      plot(tr_2,'branchlabels',true) 
%
%      % The mid-point of the original tree was the root, since it was an
%      % ultrametric tree
%      tr_3 = reroot(tr_2)
%      plot(tr_3,'branchlabels',true) 
%       
%   See also PHYTREE, PHYTREE/GET, PHYTREE/GETBYNAME, PHYTREE/PRUNE,
%   PHYTREE/SELECT, SEQNEIGHJOIN. 

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.8.3 $ $Author: batserve $ $Date: 2006/06/16 20:06:46 $

if numel(tr)~=1
     error('Bioinfo:phytree:reroot:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end
numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 

if nargin == 1;
    [node,distance] = midpoint(tr);
else
    % validate node
    if islogical(node)
        if any(numel(node) == [numLabels numLeaves])
            node = find(node);
        elseif numel(node) == numBranches
            node = find(node) + numLeaves;
        else
            error('Bioinfo:phytree:reroot:IncorrectSizeInputVector',...
                'Logical vector must have the same number of elements as nodes in the Phylogenetic Tree.');
        end
    end
    if ~isscalar(node) || node<1 || node>numLabels
        error('Bioinfo:phytree:reroot:InvalidInputNode',...
            'Invalid value for NODE.');
    end
    node = round(node);
    % when no distance is given put the root at half the branch
    if nargin<3
        distance = tr.dist(node)/2;
    end
end

% find parents for every tree node
parent = zeros(numLabels,1);
parent(tr.tree) = repmat(numLeaves+1:numLabels,1,2);

% validate distance, if necessary shift the origin node
if ~isscalar(distance) || distance<0
    error('Bioinfo:phytree:reroot:InvalidInputDistance',...
        'Invalid value for DISTANCE.');
end

% validate distance, if necessary shift the origin node
while (tr.dist(node)<=distance) && (node ~= numLabels)
    distance = distance - tr.dist(node);
    node = parent(node);
end

if node == numLabels 
    if distance>0
        warning('Bioinfo:phytree:reroot:BeyondRoot',...
                'Distance goes beyond the root, tree unchanged.')
    else
        tr.names{end} = 'Root';
    end
    return
end

% check if we just need to move the branches of current root
if parent(node) == numLabels
    tr.dist(setxor(tr.tree(end,:),node)) = ...
                sum(tr.dist(tr.tree(end,:))) - distance;
    tr.dist(node) = distance;
    tr.names{end} = 'Root';
    return
end
    
% path to root from node and bros for every point in the path
path2r = false(numBranches,1);
pathBros   = [];
me = node;
par = parent(node);
while par
    pathBros = [pathBros;setxor(tr.tree(par-numLeaves,:),me)];
    path2r(par-numLeaves) = true;
    me = par;
    par = parent(par);
end
path2rInd=find(path2r)+numLeaves;

% new tree pointers
tr.tree = [tr.tree(~path2r,:);...
      [[pathBros(end-1:-1:1);node],[pathBros(end);path2rInd(end-1:-1:1)]]];
       
% swapping distances in the nodes that belong to the path2r           
tr.dist(pathBros(end)) = tr.dist(pathBros(end)) + tr.dist(path2rInd(end-1));
tr.dist(path2rInd(2:end-1)) = tr.dist(path2rInd(1:end-2));
tr.dist(parent(node)) = tr.dist(node) - distance;
tr.dist(node) = distance;

% some branches changed positions, need to tree 
permuta = [(1:numLeaves)';find(~path2r)+numLeaves;...
            path2rInd(end-1:-1:1);numLabels];
ipermuta(permuta)=1:numLabels;
tr.tree = ipermuta(tr.tree);
tr.dist = tr.dist(permuta);
tr.names = tr.names(permuta);
tr.names{end} = 'Root';

% re-order leaves for no branch crossings
tr = prettyorder(tr);

%-% ----------------------------------------------------------------
% Selects the point where the mean of the branch length is equalized
function [node,distance] = midpoint(tr)

numBranches = size(tr.tree,1);
numLeaves = numBranches + 1;
numLabels = numBranches + numLeaves; 


branchWidth = ones(numLabels,1);
downDist = zeros(numLabels,1);
upDist = zeros(numLabels,1);
% cumulative distance and width downwards the tree
for ind = 1:numBranches
    branchWidth(numLeaves+ind) = sum(branchWidth(tr.tree(ind,:)));
    downDist(numLeaves+ind) = sum(downDist(tr.tree(ind,:)) + ...
              tr.dist(tr.tree(ind,:)).*branchWidth(tr.tree(ind,:)));
end 
% backpropagate distances 
for ind = numBranches:-1:1
    upDist(tr.tree(ind,:)) = downDist(tr.tree(ind,[2 1])) +...
       upDist(ind+numLeaves) + ...
       tr.dist(ind+numLeaves).*(numLeaves-branchWidth(ind+numLeaves)) + ...
       tr.dist(tr.tree(ind,[2 1])).*(branchWidth(tr.tree(ind,[2 1])));
end

% find all possible midpoints, solve this eq for every edge
%  ud/Nu + (x)*e = dd/Nd + (1-x)*e
%     ud = cumulative upwards distances
%     dd = cumulative downwards distances
%     Nu = number of leaves in the upper braches
%     Nd = number of leaves in the lower braches
%     e = distance of current edge

h = tr.dist~=0; % root can not be in an edge which length is zero
h(numLabels) = false; % the route can not be segmented
x = inf(numLabels,1);
x(h) = (upDist(h)./branchWidth(h) - ...
        downDist(h)./(numLeaves-branchWidth(h)))./tr.dist(h)/2 + 1/2;
    
x(h) = (upDist(h)./(numLeaves-branchWidth(h)) - ...
        downDist(h)./(branchWidth(h)))./tr.dist(h)/2 + 1/2;

% find all possible roots
h = find(x>=0 & x<=1);
if isempty(h)
    [dummy,h] = min(abs(x-1/2)); %#ok
end
% pick the most balanced one
[d,g]=min(abs(branchWidth(h)*2-numLeaves));   %#ok

node = h(g);
ratio = min(max(x(h(g)),0),1);

% if ratio is 1 then better pick the parent
if ratio == 1
   [node,dummy] = find(tr.tree==node); %#ok
   node = node + numLeaves;
   ratio = 0;
end

% change the ratio (x) to the distance to the selected node
distance = ratio .* tr.dist(node);
    
