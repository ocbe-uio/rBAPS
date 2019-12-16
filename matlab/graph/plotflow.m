function plotflow(action)
switch action
    case 'rename'
        rename;
    case 'prune'
        prune;
    case 'edit_pop_name'
        changePopNumber(0);
    case 'cancel_pop_name'
        closereq;
    case 'load_pop_names'
        loadPopNames;
    case 'ok_pop_name'
        okPopName;
    case 'five_back'
        changePopNumber(-5);
    case 'one_back'
        changePopNumber(-1);
    case 'one_ahead'
        changePopNumber(1);
    case 'five_ahead'
        changePopNumber(5);
    case 'help'
        openHelp;
end
return

% -------------------------------------------------------------------------
function rename
h0 = gcf;
h1 = findobj(h0, 'Tag','attr_menu');
g  = get(h1,'Userdata');
varnames = g.varnames;
openInputPopNamesFigure(varnames); % rename clusters

% -------------------------------------------------------------------------
function prune
h0 = gcf;
h1 = findobj(h0, 'Tag','attr_menu');
g  = get(h1,'Userdata');
adjmat = g.adjmat;
if isfield(g,'adjmat2')
    adjmat2 = g.adjmat2;
    min_strength = min(adjmat2(logical(adjmat2>0)));
else
    min_strength = min(adjmat(logical(adjmat>0)));
end

%waitALittle;
answer = inputdlg( ['Specify the minimal strength of the gene flow graph such that'...
                ' edges with lower weight will be pruned'],...
    'Prune the graph',1, {num2str(min_strength)});
if isempty(answer)  % cancel has been pressed
    return
else
    min_strength = str2num(answer{1});
    fprintf('Minimal flow strength: %4.5f\n', min_strength);
    adjmat(adjmat<min_strength) = 0;
    g.adjmat2 = adjmat;
    set(h1,'Userdata',g);
    
    % redraw the figure
    handle = gcf;
    % close(handel)
    graphvizpath = g.graphvizpath;
    nodecolors = g.nodecolors;
    nodestyles = g.nodestyles;
    arccolors = g.arccolors;
    arcstyles = g.arcstyles;
    groupnames = g.varnames;

    d = cd;
    plotmodel(g.adjmat2,g.k,'graphvizpath', graphvizpath, ...
        'nodecolors', nodecolors, ...
        'nodestyles', nodestyles, ...
        'arccolors',arccolors, ...
        'arcstyles',arcstyles, ...
        'varnames', groupnames, ...
        'target', ['matlab:',num2str(handle)]);
    % drawnow;
    cd(d);
end


%--------population_names_figure's CALLBACKS:------------------------------

function openInputPopNamesFigure(popname)
%Opens popultion_names_figure:
h1 = population_names_figure;
% h0 = findobj('Tag','input_pop_name_text');
% set(h0,'String','Name of Cluster 1:');

%Set the name of the first population to the screen:
h0 = findobj('Tag','pop_name_edit');
set(h0,'String',popname{1});

clear tiedot;
tiedot.tempnamesUD = popname;
tiedot.currentnumberUD = 1;
set(h1,'UserData',tiedot);

disableMovingButtonsIfNeeded(1);

function changePopNumber(num_of_moves)
%Updates the text, edit-field and UserData of
%population_names_figure, when some of the moving buttons has
%been pushed. Also disables some moving buttons if needed.
%Saves the name of the population of the previous number.
%This function shouldn't be called if the movement specified by
%num_of_moves isn't possible.

%Find out information of the current situation:
h1 = findobj('Tag','population_names_figure');
tiedot = get(h1,'UserData');
popname = tiedot.tempnamesUD;
current = tiedot.currentnumberUD;
%Get the old name from the 'pop_name_edit':
h0 = findobj('Tag','pop_name_edit');
name = get(h0,'String');
popname{current} = name;
%Update current and set text and edit-field accordingly:
current = current + num_of_moves;
h0 = findobj('Tag','input_pop_name_text');
set(h0,'String',['Name of Cluster ' num2str(current) ':']);
h0 = findobj('Tag','pop_name_edit');
set(h0,'String',popname{current});
%Disable movingbuttons if needed:
disableMovingButtonsIfNeeded(current);
%Save new information to 'population_names_figure':s UserData:
tiedot.tempnamesUD = popname;
tiedot.currentnumberUD = current;
set(h1,'UserData',tiedot);


function disableMovingButtonsIfNeeded(num_of_curr)
%i.e. if there are 12 populations, and num_of_curr == 9,
%five_ahead must be disabled because it takes out of the range.

%Find out the number of populations:
h0 = findobj('Tag','population_names_figure');
tiedot = get(h0,'UserData');
pops = tiedot.tempnamesUD;
num_of_pops = length(pops);
%Disable buttons if needed, enable otherwise:
if (num_of_curr - 5) < 1
    disable('five_back_button');
else enable('five_back_button');
end;
if (num_of_curr - 1) < 1
    disable('one_back_button');
else enable('one_back_button');
end;
if (num_of_curr + 5) > num_of_pops
    disable('five_ahead_button');
else enable('five_ahead_button');
end;
if (num_of_curr + 1) > num_of_pops
    disable('one_ahead_button');
else enable('one_ahead_button');
end;


function loadPopNames
%Loads the population names to UserData of 'population_
%names_figure'. Gets the names from a file specified by
%the user.

[filename,pathname] = uigetfile('*.txt','Load Population Names');
if (filename == 0) & (pathname == 0)
    %Cancel was pushed.
    return,
end;
input_file = [pathname filename];
%Read population names from the file to a variable 'names':
fid = fopen(input_file);
if fid == -1
    %File didn't exist
    msgbox('Loading of the population names was unsuccessful', ...
        'Error', 'error');
    return;
end;
line = fgetl(fid);
counter = 1;
while (line ~= -1) && ~isempty(line)
    names{counter} = line;
    line = fgetl(fid);
    counter = counter + 1;
end;
fclose(fid);
%Check that the number of lines is the same as npops:
h0 = findobj('Tag','data_info_button');
tiedot = get(h0,'UserData');
counts = tiedot.countsUD;
sizc = size(counts); npops = sizc(3);
length_names = length(names);
if npops ~= length_names
    msgbox(['Loading of names was unsuccessful.' ...
        'The number of lines in a file that contains the names ' ...
        'must be same as the number of observed ' ...
        'sampling units in the data.'] ,'Error', ...
        'error');
    return;
end;
%Save names to the UserData of the 'population_names_figure':
h0 = findobj('Tag','population_names_figure');
tiedot = get(h0,'UserData');
tiedot.tempnamesUD = names;
set(h0,'UserData',tiedot);
%Update text-field to initial state:
h0 = findobj('Tag','input_pop_name_text');
set(h0,'String','Name of Cluster 1:');
%Set the name of the first population to the screen:
h0 = findobj('Tag','pop_name_edit');
set(h0,'String',names{1});


function okPopName
%Saves the names of populations and closes
%'population_names_figure'.

h0 = findobj('Tag','population_names_figure');
tiedot = get(h0,'UserData');
names = tiedot.tempnamesUD;
%Close the figure:
% closereq;
close(h0);

handle = gcf;
h0 = findobj(handle,'Tag','attr_menu');
g  = get(h0,'Userdata');
g.varnames = names;
set(h0,'UserData',g);

% redraw the figure
switch g.type
    case 'GENEFLOW'
        graphvizpath = g.graphvizpath;
        nodecolors = g.nodecolors;
        nodestyles = g.nodestyles;
        arccolors = g.arccolors;
        arcstyles = g.arcstyles;
        groupnames = g.varnames;

        d = cd;
        plotmodel(g.adjmat2,g.k,'graphvizpath', graphvizpath, ...
            'nodecolors', nodecolors, ...
            'nodestyles', nodestyles, ...
            'arccolors',arccolors, ...
            'arcstyles',arcstyles, ...
            'varnames', groupnames, ...
            'target', ['matlab:',num2str(handle)]);
        cd(d);
    case {'nj', 'upgma'}
        switch g.visualtype
            case {'square', 'angular'}
                viewDendrogram(g.visualtype)
            case {'radial', 'phylogram'}
                viewUnrooted(g.visualtype)
        end
    otherwise
        return
end

function enable(obj_tag)
h0 = findobj('Tag',obj_tag);
set(h0,'Enable','on');


function disable(obj_tag)
h0 = findobj('Tag',obj_tag);
set(h0,'Enable','off');


function openHelp
info{1}='Gene Flow Between Populations';
info{2}='';
info{3}='Produced using GraphViz';
info{4}='Source: gene flow matrix produced by BAPS';
info{5}='';
info{6}='Author: Jing Tang';
info{7}='';
helpdlg(info,'Help');
