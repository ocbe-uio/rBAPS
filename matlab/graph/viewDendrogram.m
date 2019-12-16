function viewDendrogram(action)
% VIEWDENDROGRAM function to called by VIEWNJ to draw dendrogram
%  trees

handle = gcf;
h0 = findobj(handle,'Tag','attr_menu');
g  = get(h0,'Userdata');
g.visualtype = action;
set(h0,'Userdata',g); % store in the attribute menu
cla
axis off
% clf
% close(g.handle);
% h0 = figure('NumberTitle','off');
% g.handle = h0;
% set(h0,'menubar','none','toolbar','figure');
% set(h0,'Tag','nj_plot');
% set(h0,'Name',['Neighbor-Joining tree - ' g.filename]);;
if strcmp(g.type,'NJ')
    t = seqNeighJoin(g.D, 'equivar', g.varnames);
else 
    t = seqlinkage(g.D, 'average', g.varnames);
end

switch action
    case 'square'
       % plotNJ(g.D, 1, char(g.varnames));       
       Plot(t,'type','square');
    case 'angular'
       % plotNJ(g.D, 0, char(g.varnames)); 
       Plot(t,'type','angular');
end

% h1 = uimenu('Parent',handle, ...
%     'Label','Attributes', ...
%     'Tag','attr_menu');
% h2 = uimenu('Parent',h1, ...
%     'Label','Rename clusters', ...
%     'callback', 'plotflow rename', ...
%     'Tag','clustername_menu');
% h2 = uimenu('Parent',h1, ...
%     'Label','Visual type', ...
%     'Tag','visualtype_menu');
% h3 = uimenu('Parent',h2, ...
%     'Label','Square', ...
%     'callback', 'viewDendrogram(''square'')', ...
%     'Tag','viewsquare_menu');
% h3 = uimenu('Parent',h2, ...
%     'Label','Angular', ...
%     'callback', 'viewDendrogram(''angular'')', ...
%     'Tag','viewangular_menu');
% h3 = uimenu('Parent',h2, ...
%     'Label','Radial', ...
%     'callback', 'viewUnrooted', ...
%     'Tag','viewradial_menu');
%set(h1,'Userdata',g); % store in the attribute menu