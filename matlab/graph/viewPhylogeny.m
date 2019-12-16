function viewPhylogeny(action)
% VIEWPHYLOGENY function draws phylogenetic trees

% load the mixture result

warning('off','MATLAB:dispatcher:InexactMatch')

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
h0 = findobj('Tag','filename1_text');
filename = get(h0,'String');

aln.npops = c.npops;
aln.COUNTS = c.COUNTS;
aln.adjprior = c.adjprior;
aln.data = c.data;
aln.partition = c.PARTITION;
aln.clusternames = ''; % default cluster names are 1 x npops integers.
aln.noalle = c.noalle;
clear c;

disp('---------------------------------------------------');
if strcmp(action, 'nj')
    disp('Viewing the Neighbor-Joining tree.');
else
    disp('Viewing the UPGMA tree.');
end
disp(['Load the mixture result from: ',[filename],'...']);

% Select the type of distance matrix
[D] = chooseDistance(aln);
if (isempty(D))
    return;
else
    % draw the window
    h0 = figure('NumberTitle','off');
    set(h0, 'NumberTitle', 'off');
    set(h0,'menubar','none','toolbar','figure');
    set(h0,'Tag','nj_plot');
    set(h0,'Name',[ upper(action) ' tree - ' filename]);
    h1 = uimenu('Parent',h0, ...
        'Label','Attributes', ...
        'Tag','attr_menu');
    h2 = uimenu('Parent',h1, ...
        'Label','Rename clusters', ...
        'callback', 'plotflow rename', ...
        'Tag','clustername_menu');
    h2 = uimenu('Parent',h1, ...
        'Label','Visual type', ...
        'Tag','visualtype_menu');
    h3 = uimenu('Parent',h2, ...
        'Label','Square', ...
        'callback', 'viewDendrogram(''square'')', ...
        'Tag','viewsquare_menu');
    h3 = uimenu('Parent',h2, ...
        'Label','Angular', ...
        'callback', 'viewDendrogram(''angular'')', ...
        'Tag','viewangular_menu');
    h3 = uimenu('Parent',h2, ...
        'Label','Radial', ...
        'callback', 'viewUnrooted(''radial'')', ...
        'Tag','viewradial_menu');
    h3 = uimenu('Parent',h2, ...
        'Label','Phylogram', ...
        'callback', 'viewUnrooted(''phylogram'')', ...
        'Tag','viewphylogram_menu');

    if strcmp(action, 'nj')
        t = seqNeighJoin(D, 'equivar', correct(aln.clusternames, aln.npops));
    else
        t = seqlinkage(D, 'average', correct(aln.clusternames, aln.npops));
    end
    Plot(t);
end

% save the parameters
% g.handle = h0;
g.varnames = correct(aln.clusternames, aln.npops);
g.D = D;
g.type = action;
g.filename = filename;
g.visualtype = 'square';
g.tree = t;
set(h1,'Userdata',g); % store in the attribute menu

%--------------------------------------------------------------------------
% SUBFUNCTIONS
%--------------------------------------------------------------------------
function varnames = correct(clusternames,npops)

varnames = cell(1,npops);
if isempty(clusternames)
    for i=1:npops
        varnames{i} = sprintf('Cluster %d',i);
    end
end
