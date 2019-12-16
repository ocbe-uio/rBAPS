function plotPhytree(tr,varargin)
%PLOTPHYTREE renders a phylogenetic tree.
%
%   PLOTPHYTREE(TREE) renders a phylogenetic tree object into a MATLAB figure as a
%   phylogram. 
%
%   plotphytree(...,'ROTATION',value) will orient the phylogenetic tree within
%   the figure window. Positive angles cause counterclockwise rotation,
%   otherwise clockwise rotation.   
%
%   plotphytree(...,'FONTSIZE',value) will set the lable font size.A value
%   specifying the font size to use for text in units determined by the FontUnits
%   property (1 point = 1/72 inch). The default is calculated from the data
%   according to the number of nodes.
%
%   plotphytree(...,'LINESTYLE',value) will set the color of the tree. The
%   default color is blue.
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot 
%          c     cyan          +     plus               --    dashed   
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%                              v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram
%
%   plotphytree(...,'FONTCOLOR',value) will set the color of lable. It is a
%   vector includes [R, G, B] components. Each with the range of 0 to 1.
%   The default setting is [.2 .2 .2].
%
%   
%   Example:
%      
%       tr = phytreeread('pf00002.tree');
%       plotphytree(tr,'ROTATION',-pi/2, 'FONTSIZE', 8, 'LINESTYLE', 'b', 'FONTCOLOR', [0.2 0.2 0.2]);
%
%       

if numel(tr)~=1
     error('plotphytree:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

% set defaults
rotation = 0;
fontsize = 1;
fontsizeset = false;
linestyle = '-b';
fontcolor = [.2 .2 .2];


if nargin>1 && islogical(varargin{1})
    activeBranches = varargin{1};  
    argStart = 2;
else     
    argStart = 1;
end

if nargin - argStart > 0
        if rem(nargin - argStart,2) == 1
        error('plotphytree:IncorrectNumberOfArguments',...
              'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'rotation', 'fontsize', 'linestyle', 'fontcolor'}; 
    
    for j = argStart:2:nargin-argStart          
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('plotphytree:UnknownParameterName',...
                  'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('plotphytree:AmbiguousParameterName',...
                  'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % rotation
                    if isreal(pval(1))
                        rotation = double(pval(1));
                    else
                        error('plotphytree:NotValidType',...
                              'ROTATION must be numeric and real');
                    end
                case 2 % fontsize
                    if isreal(pval(1))
                        fontsize = uint8(pval(1));
                        fontsizeset = true;
                    else
                        error('plotphytree:NotValidType',...
                              'fontsize must be numeric');
                    end  
                case 3 % linecolor
                    linestyle = pval(1);
                case 4 % Fontcolor
                    fontcolor = [pval(1) pval(2) pval(3)];                    
            end
        end
    end
end    

orgtr=tr;
treeinfo = GET(tr);

%Get tree structure
tr = struct(tr);
tr.numLeaves = treeinfo.NumLeaves;
tr.numNodes = treeinfo.NumNodes;
tr.numBranches = treeinfo.NumBranches;
tr.LeafNames = treeinfo.LeafNames;
tr.BranchNames = treeinfo.BranchNames;
tr.numTree = size(tr.tree, 1);

% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% obtain numbers of leaves for each nodes
tr.numChildrenofNode = ChildCount(tr.par, tr.numLeaves, tr.numNodes);
tr.numChildrenofNode(tr.numChildrenofNode == 0) = 1;

% find angle of each node 
unitDegree = 2 * pi / tr.numLeaves;

tr.NodeSector = tr.numChildrenofNode * unitDegree;    
tr.NodeSector(tr.NodeSector == 0) = unitDegree;
tr.NodeSector(tr.numNodes) = 0;

tr.NodeAngle = zeros(tr.numNodes, 1);
tr.NodeAngle(tr.tree(tr.numTree, 2)) = tr.NodeSector(tr.tree(tr.numTree, 2)) / 2;
tr.NodeAngle(tr.tree(tr.numTree, 1)) = tr.NodeSector(tr.tree(tr.numTree, 2)) + tr.NodeSector(tr.tree(tr.numTree, 1)) / 2;
for i = tr.numTree - 1 : -1 :1
    tr.NodeAngle(tr.tree(i, 2)) = tr.NodeAngle(tr.par(tr.tree(i, 2))) - tr.NodeSector(tr.par(tr.tree(i, 2))) / 2 + tr.NodeSector(tr.tree(i, 2)) / 2;
    tr.NodeAngle(tr.tree(i, 1)) = tr.NodeAngle(tr.par(tr.tree(i, 1))) + tr.NodeSector(tr.par(tr.tree(i, 1))) / 2 - tr.NodeSector(tr.tree(i, 1)) / 2;
end

% Rotation
tr.NodeAngle = tr.NodeAngle + rotation;

% find (x, y) coordinates of nodes
tr.d = tr.dist; 
tr.x = tr.d .* cos(tr.NodeAngle);
tr.y = tr.d .* sin(tr.NodeAngle);
 
for i = tr.numNodes - 1: -1 : 1
    tr.x(i) = tr.x(i) + tr.x(tr.par(i));
    tr.y(i) = tr.y(i) + tr.y(tr.par(i));
end

nodeIndex = 1 : tr.numNodes;
X = tr.x([nodeIndex;[tr.par(1:tr.numNodes-1) tr.numNodes]]);
Y = tr.y([nodeIndex;[tr.par(1:tr.numNodes-1) tr.numNodes]]); 
tr.txtAngle = cal_textAngle(tr);

fig = gcf;
% fig = figure('Renderer','ZBuffer');
h.fig = fig;
h.axes = axes; hold on;
sepUnit = max(tr.x)*[-1/20 21/20];

set(h.axes,'XTick',[],'YTick',[]);
set(h.axes,'Position',[.05 .05 .9 .9])
dispTerminalLabels = false;
axis equal

h.BranchLines = plot(X,Y,linestyle);

% resize figure if needed
temp = 10/pi*tr.numLeaves;
correctFigureSize(fig,temp,temp);
fontRatio = max(get(fig,'Position').*[0 0 1 0])/tr.numLeaves;
set(h.axes, 'Fontsize', fontsize);

% set leaf nodes labels
leafIndex = 1 : tr.numLeaves;
X = tr.x(leafIndex);
Y = tr.y(leafIndex);

txtangle = tr.txtAngle;

% % Rotate Label
h.leafNodeLabels = zeros(1, tr.numLeaves);
for i = 1:tr.numLeaves
    h.leafNodeLabels(i) = text(X(i),Y(i),tr.names(i), 'Rotation', txtangle(i));
end

set(h.leafNodeLabels,'color',fontcolor,'clipping','on')
if fontsizeset
    set(h.leafNodeLabels, 'Fontsize', fontsize);
else
    set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));
end
 
textHeight = mean(cell2mat(get(h.leafNodeLabels,'Extent')))*[0 0 0 1]';
for ind = 1:numel(h.leafNodeLabels)
    if X(ind) - tr.x(tr.par(ind)) < 0
        if txtangle(ind) < 90.0 - 0.001 || txtangle(ind) > 90.0 + 0.001
            set(h.leafNodeLabels(ind),'horizontal','right')
            set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[sepUnit(1)*cos(txtangle(ind) * 2 * pi/360) sepUnit(1)*sin(txtangle(ind) * 2 * pi/360) 0])
        else
            if Y(ind) > 0 
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[0 sepUnit(1) 0])
            else
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[0 sepUnit(1) 0])
            end
        end
    else
        if txtangle(ind) < 90.0 - 0.001 || txtangle(ind) > 90.0 + 0.001
            set(h.leafNodeLabels(ind),'horizontal','left')
            set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[sepUnit(1)*cos(txtangle(ind) * 2 * pi/360) sepUnit(1)*sin(txtangle(ind) * 2 * pi/360) 0])
        else
            if Y(ind) > 0 
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[0 sepUnit(1) 0])
            else
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[0 sepUnit(1) 0])
            end
        end
    end
end

dispLeafLabels = 1;
% correct axis limits given the extent of labels
if dispLeafLabels
    E = cell2mat(get(h.leafNodeLabels,'Extent'));
    if strcmp(get(gca,'XDir'),'reverse')
        E(:,1) = E(:,1) - E(:,3);
    end
    if strcmp(get(gca,'YDir'),'reverse')
        E(:,2) = E(:,2) - E(:,4);
    end
    E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
    mins = min(E(:,[1 2]));
    maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
    axis([mins(1) maxs(1) mins(2) maxs(2)])
end

if dispTerminalLabels
    set(h.terminalNodeLabels,'Fontsize',min(9,ceil(fontRatio/1.5)));  
end

box off
hold off

% store handles
set(fig,'UserData',h)   
if nargout
    handles = h;
end


%**************************************************************
% Count number of children for each node
function numChildren = ChildCount(par, numleafs, numnodes)
n = length(par);
numChildren = zeros(1, numnodes);
NoBranchStart = numleafs + 1;

for i = 1 : numleafs
    numChildren(par(i)) = numChildren(par(i)) + 1;
end

for i = NoBranchStart : n
    numChildren(par(i)) = numChildren(par(i)) + numChildren(i);
end
%**************************************************************

%*********************************************************************
function correctFigureSize(fig,recommendedHeight,recommendedWidth)
% helper function to increase initial figure size depending on the screen &
% tree sizes
screenSize = diff(reshape(get(0,'ScreenSize'),2,2),[],2)-[0;100];
            % 100 gives extra space for the figure header and win toolbar
position = get(fig,'Position');
if recommendedHeight > position(4)
    if recommendedHeight < sum(position([2 4]))
        position(2) = sum(position([2 4])) - recommendedHeight;
        position(4) = recommendedHeight;
    elseif recommendedHeight < screenSize(2)
        position(2) = 30; 
        position(4) = recommendedHeight;
    else 
        position(2) = 30; 
        position(4) = screenSize(2);
    end
end
if recommendedWidth > position(3)
    if recommendedWidth < sum(position([1 3]))
        position(1) = sum(position([1 3])) - recommendedWidth;
        position(3) = recommendedWidth;
    elseif recommendedWidth < screenSize(1)
        position(1) = 0; 
        position(3) = recommendedHeight;
    else 
        position(1) = 0; 
        position(3) = screenSize(1);
    end
end    
set(fig,'Position',position)
%*********************************************************************

%*********************************************************************
function [lefttree, righttree, thirdtree] = findsubtreeofnode(tr)

lefttree = zeros(tr.numBranches, tr.numNodes);
righttree = zeros(tr.numBranches, tr.numNodes);
thirdtree = zeros(tr.numBranches, tr.numNodes);

% thirdtree_left = zeros(tr.numBranches, tr.numNodes);
% thirdtree_right = zeros(tr.numBranches, tr.numNodes);

lefttree(1, 1) = tr.tree(1, 1);
righttree(1,1) = tr.tree(1, 2);
lefttree(1, tr.numNodes) = 1;
righttree(1, tr.numNodes) = 1;

for i = 2 : tr.numBranches
    j = 0; 
    tmp = tr.tree(i, 1);
    lefttree(i, 1) = tmp;
    lefttree(i, tr.numNodes) = 1;
    if tmp > tr.numLeaves
        lefttrlen = lefttree(tmp - tr.numLeaves, tr.numNodes);
        sub = lefttree(tmp - tr.numLeaves, 1 : lefttrlen);
        lefttree(i, 2 : lefttrlen + 2 - 1) = sub;
        
        righttrlen = righttree(tmp - tr.numLeaves, tr.numNodes);
        sub = righttree(tmp - tr.numLeaves, 1 : righttrlen);
        lefttree(i, lefttrlen + 3 - 1: lefttrlen + 3 - 1 + righttrlen - 1) = sub;  
        lefttree(i, tr.numNodes) = lefttrlen + righttrlen + lefttree(i, tr.numNodes);
    end
        
    j = 0;
    tmp = tr.tree(i, 2);
    righttree(i, 1) = tmp;
    righttree(i, tr.numNodes) = 1; 
    if tmp > tr.numLeaves
        lefttrlen = lefttree(tmp - tr.numLeaves, tr.numNodes);
        sub = lefttree(tmp - tr.numLeaves, 1 : lefttrlen);
        righttree(i, 2 : lefttrlen + 2 - 1) = sub;
        
        righttrlen = righttree(tmp - tr.numLeaves, tr.numNodes);
        sub = righttree(tmp - tr.numLeaves, 1 : righttrlen);
        righttree(i, lefttrlen + 3 - 1: lefttrlen + 3 - 1 + righttrlen - 1) = sub;  
        righttree(i, tr.numNodes) = lefttrlen + righttrlen + righttree(i, tr.numNodes);
    end        
end

% Third tree
node = 1 : tr.numLeaves;
for i = 1 : tr.numBranches
    tmp = [lefttree(i, 1:lefttree(i, tr.numNodes)), righttree(i, 1:righttree(i, tr.numNodes))];
    tf = ismember(node, tmp);
    tmp = node(find(tf == 0));
    tmp = tmp(find(tmp <= tr.numLeaves));
    nlength = length(tmp);
    thirdtree(i, 1 : nlength) = tmp;
    thirdtree(i, tr.numNodes) = nlength;
end

%*********************************************************************
function [treeangle, leftnode, rightnode] = caltreeangle(treenode, root, tr)
tmp = find(treenode <= tr.numLeaves);
leftnode = treenode(tmp(1));
rightnode = treenode(tmp(length(tmp)));

x0 = tr.x(root);
y0 = tr.y(root);
x1 = tr.x(leftnode);
y1 = tr.y(leftnode);
x2 = tr.x(rightnode);
y2 = tr.y(rightnode);

treeangle = cal_angle(x0, y0, x1, y1, x2, y2);
 
%*********************************************************************
function [subtree1, subtree2] = findsubtree(node, tr)
subtree1 = [];
subtree2 = [];

nodeindex = 1 : tr.numNodes;

i1 = 0;
i2 = 0;

flag1 = 0;
flag2 = 0;    

while ~flag1 || ~flag2
    nodetmp = find(tr.par == node);
    
    if nodetmp(1) <= tr.numLeaves
        subtree1 = [subtree1, nodetmp(1)];        
        i1 = i1 + 1;
        flag1 = 1;
    else
        if i1 == 0
            subtree1 = [subtree1, nodetmp(1)];
        end
        [subtr1, subtr2] = findsubtree(nodetmp(1), tr);
        flag1 = 1;
        subtree1 = [subtree1, subtr1, subtr2];        
    end
    
    if nodetmp(2) <= tr.numLeaves
        subtree2 = [subtree2, nodetmp(2)];
        i2 = i2 + 1;
        flag2 = 1;
    else
        if i2 == 0
            subtree2 = [subtree2, nodetmp(2)];
        end
        [subtr1, subtr2] = findsubtree(nodetmp(2), tr);
        flag2 = 1;
        subtree2 = [subtree2, subtr1, subtr2];
    end
end

%*********************************************************************
function leavenodes = findanothersubtree(node, tr)
nextnode = node;
cutpoint = tr.par(length(tr.par));
subtree = [];
tmpnode = nextnode;
nextnode = tr.par(nextnode);

while nextnode <= cutpoint
    [subtree1, subtree2] = findsubtree(nextnode, tr);
    if sum(find(subtree1 == node)) == 0
        subtree = [subtree, subtree1];
    else
        subtree = [subtree, subtree2];
    end    
    tmpnode = nextnode;
    if nextnode < cutpoint
        nextnode = tr.par(nextnode); 
    else 
        break;
    end
end
leavenodes = subtree(find(subtree <= tr.numLeaves));
    
%*********************************************************************
function slope = cal_slope(x1, y1, x2, y2)
if x1 == x2
    slope = tan(pi/2);
else
    slope = (y2 - y1)/(x2 - x1);
end

%*********************************************************************
function xangle = cal_angle(x0, y0, x1, y1, x2, y2)
k1 = cal_slope(x0, y0, x1, y1);
k2 = cal_slope(x0, y0, x2, y2);
xangle = atan(abs((k2 - k1)/(1 + k1 * k2)));

%*********************************************************************
function [leftnode, rightnode, xangle] = cal_thr_angle(trnode, node, tr)
leftnode = 0;
rightnode = 0; 
xangle = 0;

n = length(trnode);
if n == 1
    leftnode = trnode(1);
    rightnode = trnode(1);
    xangle = 0;
    return;
end

x0 = tr.x(node);
y0 = tr.y(node);

for i = 1 : n - 1
    x1 = tr.x(trnode(i));
    y1 = tr.y(trnode(i));    
    k1 = cal_slope(x0, y0, x1, y1);
    for j = i + 1 : n
        x2 = tr.x(trnode(j));
        y2 = tr.y(trnode(j));
        k2 = cal_slope(x0, y0, x2, y2);
        theta = atan(abs((k2 - k1)/(1 + k1 * k2)));
        if xangle < theta
            xangle = theta;
            leftnode = trnode(i);
            rightnode = trnode(j);
        end
    end
end

%*********************************************************************
function [left, right] = leftnrightnode(subtree, tr)
lefttree = tr.tree(:, 1);
righttree = tr.tree(:, 2);

left = 0;
right = 0;
numindex = length(lefttree);
i = 1;
for i = 1 : numindex
    [tf, index] = ismember(lefttree(i), subtree);
    if tf == 1 
        left = lefttree(i);
        break;
    end
end

for i = 1 : numindex
    [tf, index] = ismember(righttree(i), subtree);
    if tf == 1
        right = righttree(i);
        break;
    end
end

if left == 0
    left = right;
end
if right == 0 
    right = left;
end

%*********************************************************************
function nodeangle = cal_nodeangle(y, x)
nodeangle = atan2(y, x);
if nodeangle < 0
    nodeangle = 2*pi + nodeangle;
end

%*********************************************************************
function tree = adjustsubtreeangle(subtr1, subtr2, subtr3, orgnode, tr)

tr1node = subtr1(find(subtr1 <= tr.numLeaves));
tr2node = subtr2(find(subtr2 <= tr.numLeaves));

tr.x(subtr1) = tr.x(subtr1) - tr.x(orgnode);
tr.y(subtr1) = tr.y(subtr1) - tr.y(orgnode);

tr.x(subtr2) = tr.x(subtr2) - tr.x(orgnode);
tr.y(subtr2) = tr.y(subtr2) - tr.y(orgnode);

tr.x(subtr3) = tr.x(subtr3) - tr.x(orgnode);
tr.y(subtr3) = tr.y(subtr3) - tr.y(orgnode);

[tree1left, tree1right] = leftnrightnode(tr1node, tr);
[tree2left, tree2right] = leftnrightnode(tr2node, tr);
tr3node = subtr3;%leftnrightnode(subtr3, tr);
tree3left = tr3node(1);
tree3right = tr3node(2);    

a1left = cal_nodeangle(tr.y(tree1left), tr.x(tree1left));
a1right = cal_nodeangle(tr.y(tree1right), tr.x(tree1right));

a2left = cal_nodeangle(tr.y(tree2left), tr.x(tree2left));
a2right = cal_nodeangle(tr.y(tree2right), tr.x(tree2right));

a3left = cal_nodeangle(tr.y(tree3left), tr.x(tree3left));
a3right = cal_nodeangle(tr.y(tree3right), tr.x(tree3right));

daylight = [a1right - a3left, a2right - a1left, a3right - a2left];
equaldaylight = sum(abs(daylight))/3;

adjust1n3 = equaldaylight - daylight(1);
adjust2n3 = daylight(3) - equaldaylight;

% Rotation --- subtr1
a = tr.x(subtr1);
b = tr.y(subtr1);
if adjust1n3 ~= 0
    tr.x(subtr1) = a * cos(adjust1n3) - b * sin(adjust1n3);
    tr.y(subtr1) = a * sin(adjust1n3) + b * cos(adjust1n3);
end

% Rotation --- subtr2
a = tr.x(subtr2);
b = tr.y(subtr2);
if adjust2n3 ~= 0
    tr.x(subtr2) = a * cos(adjust2n3) - b * sin(adjust2n3);
    tr.y(subtr2) = a * sin(adjust2n3) + b * cos(adjust2n3);
end

tr.x(subtr1) = tr.x(subtr1) + tr.x(orgnode);
tr.y(subtr1) = tr.y(subtr1) + tr.y(orgnode);

tr.x(subtr2) = tr.x(subtr2) + tr.x(orgnode);
tr.y(subtr2) = tr.y(subtr2) + tr.y(orgnode);

tr.x(subtr3) = tr.x(subtr3) + tr.x(orgnode);
tr.y(subtr3) = tr.y(subtr3) + tr.y(orgnode);

tree = tr;

%*********************************************************************
function tree = adjusttree(lefttree, righttree, thirdtree, node, tr)

tr1 = lefttree(node - tr.numLeaves, 1 : lefttree(node - tr.numLeaves, tr.numNodes));
tr2 = righttree(node - tr.numLeaves, 1 : righttree(node - tr.numLeaves, tr.numNodes));
tr3 = thirdtree(node - tr.numLeaves, 1 : thirdtree(node - tr.numLeaves, tr.numNodes));
[lefttreeangle, l_leftnode, l_rightnode] = caltreeangle(tr1, node, tr);
[righttreeangle, r_leftnode, r_rightnode] = caltreeangle(tr2, node, tr);
[thr_leftnode, thr_rightnode, xangle] = cal_thr_angle(tr3, node, tr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1 = sqrt((tr.x(l_rightnode) - tr.x(thr_leftnode)) * (tr.x(l_rightnode) - tr.x(thr_leftnode)) + (tr.y(l_rightnode) -tr.y(thr_leftnode)) * (tr.y(l_rightnode) -tr.y(thr_leftnode)));
d2 = sqrt((tr.x(l_rightnode) - tr.x(thr_rightnode)) * (tr.x(l_rightnode) - tr.x(thr_rightnode)) + (tr.y(l_rightnode) -tr.y(thr_rightnode)) * (tr.y(l_rightnode) -tr.y(thr_rightnode)));

if d2 > d1
    tmp = thr_leftnode;
    thr_leftnode = thr_rightnode;
    thr_rightnode = tmp;
end

% daylight = 2 * pi - lefttreeangle - righttreeangle - xangle;
% equaldaylight = daylight / 3;

[angle1to3, angle2to3, angle1to2] = cal_daylight(l_leftnode, l_rightnode, r_leftnode, r_rightnode, thr_leftnode, thr_rightnode, node, tr);

daylight = angle1to3 + angle2to3 + angle1to2;
equaldaylight = daylight / 3;

adjust1n3 = equaldaylight - angle1to3;
adjust2n3 = angle2to3 - equaldaylight;

tr.x(tr1) = tr.x(tr1) - tr.x(node);
tr.y(tr1) = tr.y(tr1) - tr.y(node);

tr.x(tr2) = tr.x(tr2) - tr.x(node);
tr.y(tr2) = tr.y(tr2) - tr.y(node);

% Rotation --- subtr1
a = tr.x(tr1);
b = tr.y(tr1);
if adjust1n3 ~= 0
    tr.x(tr1) = a * cos(adjust1n3) - b * sin(adjust1n3);
    tr.y(tr1) = a * sin(adjust1n3) + b * cos(adjust1n3);
end

% Rotation --- subtr2
a = tr.x(tr2);
b = tr.y(tr2);
if adjust2n3 ~= 0
    tr.x(tr2) = a * cos(adjust2n3) - b * sin(adjust2n3);
    tr.y(tr2) = a * sin(adjust2n3) + b * cos(adjust2n3);
end

tr.x(tr1) = tr.x(tr1) + tr.x(node);
tr.y(tr1) = tr.y(tr1) + tr.y(node);

tr.x(tr2) = tr.x(tr2) + tr.x(node);
tr.y(tr2) = tr.y(tr2) + tr.y(node);

tree = tr;

%*********************************************************************
function [angle1to3, angle2to3, angle1to2] = cal_daylight(l_leftnode, l_rightnode, r_leftnode, r_rightnode, thr_leftnode, thr_rightnode, node, tr)

x0 = tr.x(node);
y0 = tr.y(node);

llx = tr.x(l_leftnode);
lly = tr.y(l_leftnode);
lrx = tr.x(l_rightnode);
lry = tr.y(l_rightnode);
rlx = tr.x(r_leftnode);
rly = tr.y(r_leftnode);
rrx = tr.x(r_rightnode);
rry = tr.y(r_rightnode);

tlx = tr.x(thr_leftnode);
tly = tr.y(thr_leftnode);
trx = tr.x(thr_rightnode);
tryy = tr.y(thr_rightnode);

angle1to3 = cal_angle(x0, y0, llx, lly, trx, tryy);
angle2to3 = cal_angle(x0, y0, rrx, rry, tlx, tly);
angle1to2 = cal_angle(x0, y0, lrx, lry, rlx, rly);

%*********************************************************************
function txtAngle = cal_textAngle(tr)
nodeIndex = 1 : tr.numLeaves;
X = tr.x(nodeIndex) - tr.x(tr.par(nodeIndex));
Y = tr.y(nodeIndex) - tr.y(tr.par(nodeIndex)); 

txtAngle = atan(Y./X) * 360/(2*pi);