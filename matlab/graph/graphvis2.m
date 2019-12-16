function graphvis2(adj_mat, linkage_file_name, inliers, graphviz_path, npops)
% Filename: graphvis2.m
% graphvis2(adj_mat, linkage_file_name)
%
% Description: 
% Graph visualization given the adjacency matrix. 

% Author: Jing Tang
% Modified date: 01/06/2005

% Input:
% adj_mat: The adjacency matrix.
% linkage_file_name:

% Output:
%

% Ensure adj_mat is symmetric and square.

[n,m] = size(adj_mat);
if n ~= m, error ('Adjacency matrix must be square'), end;
%if ~all(diag(adj_mat)), error('The diagonal must be nonzero'), end;

npops_in = length(inliers);
groupnames = cell(1,npops_in);
nodecolors = cell(npops_in,1);
allnode_color = giveColors(npops);
arccolors = cell(npops_in);
adjmat = adj_mat ~= 0;
for i=1:npops_in
    groupnames{i} = sprintf('Cluster %d',inliers(i));
    nodecolors{i} = num2str(rgb2hsv(allnode_color(inliers(i),:)));
    arccolors(i,adjmat(i,:)) = {nodecolors{i}};
end
nodestyles = cell(npops_in, 1);
nodestyles(:) = {'filled'};
arcstyles = cell(npops_in);
arcstyles(adjmat) = {'filled'};

handle = plotmodel(adj_mat',[1:npops_in],'graphvizpath', graphviz_path, ...
                             'nodecolors', nodecolors, ...
                             'nodestyles', nodestyles, ...
                             'arccolors',arccolors, ...
                             'arcstyles',arcstyles, ...
                             'varnames', groupnames);

% set(handle,'menubar','none','numbertitle','off','toolbar','figure');
set(handle,'numbertitle','off','toolbar','figure');
m = findall(gcf,'type','uimenu');
set(m([1:7]),'Visible', 'off');

h1 = uimenu('Parent',handle, ...
    'Label','Attributes', ...
    'Tag','attr_menu');
h2 = uimenu('Parent',h1, ...
    'Label','Rename clusters', ...
    'callback', 'plotflow rename', ...
    'Tag','clustername_menu');
h3 = uimenu('Parent',h1, ...
    'Label','Prune edges', ...
    'callback','plotflow prune', ...
    'Tag','edge_menu');
h4 = uimenu('Parent',handle, ...
    'Label','Help', ...
    'callback', 'plotflow help', ...
    'Tag','help_menu');
% h5 = uimenu('Parent',h1, ...
%     'Callback','baps4cbf about', ...
%     'Enable','on', ...
%     'Label','About', ...
%     'Tag','about_menu');
set(handle,'Name',[' Gene flow - ' linkage_file_name ]);

% save the parameters
g.handle = handle;
g.adjmat = adj_mat';
g.adjmat2 = g.adjmat;
g.k = [1:npops_in];
g.graphvizpath = graphviz_path;
g.nodecolors = nodecolors;
g.nodestyles = nodestyles;
g.arccolors = arccolors;
g.arcstyles = arcstyles;
g.varnames = groupnames;
g.type = 'GENEFLOW';
set(h1,'Userdata',g); % store in the attribute menu



% Old version 

% Decide the coordinate of each node.
% [p,p,r] = dmperm(adj_mat);
% nnodes = length(adj_mat);
% nblocks = length(r)-1;
% [B,ix]=sort(p,2);
% Coordinates=[[ix'],zeros(nnodes,1)];
% 
% % Plot the graph structure
% clf reset
% set(gcf, 'color', 'white', 'menubar', 'none', 'numbertitle','off','name', 'Graphical model')
% gplot(adj_mat, Coordinates, '-ob');
% h=findobj(gca, 'type','line');
% set(h,'markersize',10) 
% xlim([0 nnodes+1]);
% axis off
% text(Coordinates(:,1),Coordinates(:,2)-0.1, int2str(B(:)));
%drawnow
