function baps4cbf(action)

base = findobj('Tag','base_figure');
%setWindowOnTop(base,'false')

switch action

    case 'mix1_button' 
        greedyMix(-1);        

    case 'mix2_button'
        greedyPopMix;

    case 'trained_button'
        trainedMix;

    case 'spatial_button'
        spatialMixture;
        showmethebaps;
        
    case 'spatial2_button'
        spatialPopMixture;
        showmethebaps;
        
    case 'linkage_button'
        linkageMixture_speed;
        %linkageMixture_ultraspeed;
        showmethebaps;
       
    case 'admix1_button'
        admix1(-1);
        showmethebaps;
    case 'admix2_button'
        admix2;
        
    case 'compare_menu'
        compare;
    case 'compare_admix_menu'
        compare_admix;
        
    case 'load_mixture_menu'
        loadMixture;
        
    case 'load_admixture_menu'
        loadAdmixture;
        
    case 'load_spatial_menu'
        loadSpatial;
        
    case 'output_menu'
        asetaOutputFile;

    case 'remove_outputfile_menu'
        poistaOutputFile;

    case 'close_menu'
        closeFile;
        
    case 'exit_menu'
        h0 = findobj('Tag','base_figure'); delete(h0);
        h0 = findobj('Tag','image_figure'); delete(h0);

    case 'loadfigure_menu'
        loadFigure;

    case 'plot_coordinates_menu'
        plotCoordinates;

    case 'partitio_menu'
        viewPartitio;
     
    case 'admix_menu'
        viewAdmixture;
    
    case 'likelihood_menu'
        viewLoghood;
        
    case 'energy_menu'
        viewEnergy;
        
    case 'geneflow_menu'
        viewGeneflow;

    case 'voronoi_menu'
        voronoiTessellation;

    case 'varmuus_menu'
        localUncertainty;
        
    case 'changecolor_menu'
        changeColor;
        
    case 'helpdoc'
        openHelpDoc;
        
    case 'helponline'
        openHelpHtml;
        
    case 'about'
        openAboutWindow;
    
    case 'calculate_kl'
        calculateDis('KL');
    
    case 'calculate_nei'
        calculateDis('Nei');

    case 'calculate_hamming'
        calculateDis('Hamming');
        
    case 'upgma_menu'
        viewPhylogeny('upgma');
        
    case 'nj_menu'
        viewPhylogeny('nj');
        
    case 'mutationplot_menu'
        mutationPlot(-1);
        
    case 'fixk_menu'
        goToFixedK;
       
    case 'partitioncompare_menu'
        goToPartitionCompare;
end

return

%--------------------------------------------------------------------------
%KUVIEN LATAAMINEN
%--------------------------------------------------------------------------


function loadFigure
waitALittle;
[filename,pathname] = uigetfile('*.fig','Load Figure');
if (sum(filename)==0) || (sum(pathname)==0)
    return;
end
fig_file_name = [pathname filename];
open(fig_file_name);

% ----------------------------
% Old version
% ----------------------------
% % Loads previously saved figure.
% 
% waitALittle;
% [filename,pathname] = uigetfile('*.mat','Load Figure');
% if (sum(filename)==0) || (sum(pathname)==0)
%     return;
% end
% fig_file_name = [pathname filename];
% %Figure file format must be *.mat. Ensure it:
% isMat = isTheFileMatFile(fig_file_name);
% if isMat == 0
%     msgbox(['Only figures that have been saved in BAPS can be loaded in BAPS. ' ...
%         'Those figures have extension ".mat".'],'Error', ...
%         'error');
%     return;
% end;
% struct_array = load([pathname filename]);
% if isfield(struct_array,'tiedot')  %Matlab versio
%     tiedot = struct_array.tiedot;
%     if ~isfield(tiedot,'info')
%         disp('Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'info')  %Mideva versio
%     tiedot = struct_array;
% else
%     disp('Incorrect file format');
%     return;
% end
% 
% if isfield(tiedot, 'rows')
%     rows = tiedot.rows;
%     partition = tiedot.info;
%     popnames = tiedot.popnames;
%     viewPopMixPartition(partition, rows, popnames);
% else
%     popnames = tiedot.popnames;
%     info = tiedot.info;
%     if (size(info,2)>1)
%         %info on osuudet
%         osuudet = info;
%         viewPartition(osuudet,popnames);
%     else
%         info = round(info);
%         partition = info;
%         viewMixPartition(partition, popnames);
%     end
% end


function isMat = isTheFileMatFile(filename)
%Checks that the file 'filename' is of the
%*.mat format. If so, isMat = 1. Otherwise, isMat = 0.

len = length(filename);
if len < 5
    isMat = 0; return;
end;
ending = filename(end-3:end);
if isequal(ending,'.mat')
    isMat = 1;
else
    isMat = 0;
end;

%--------------------------------------------------------------------

function asetaOutputFile
waitALittle;
[filename, pathname] = uiputfile('*.txt', 'Specify output file');
if filename==0
    return
end

h0 = findobj('Tag','filename2_text');
set(h0,'String',[pathname filename]);


%---------------------------------------------------


function poistaOutputFile
h0 = findobj('Tag','filename2_text');
set(h0,'String','');


%-----------------------------------------------------

function plotCoordinates

waitALittle;
[filename, pathname] = uigetfile('*.txt', 'Load Coordinate File');
if filename==0
    return
end
X = load([pathname filename]);
if size(X,2)~=2
    disp('Incorrect file format');
    return
end

waitALittle;
[filename, pathname] = uigetfile('*.mat', 'Load mixture clustering of individuals');
%load([pathname filename],'c');
struct_array = load([pathname filename]);
if isfield(struct_array,'c')  %Matlab versio
    c = struct_array.c;
    if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd')
        disp('Incorrect file format');
        return
    end
elseif isfield(struct_array,'PARTITION')  %Mideva versio
    c = struct_array;
    if ~isfield(c,'rowsFromInd')
        disp('Incorrect file format');
        return
    end
else
    disp('Incorrect file format');
    return;
end
PARTITION = c.PARTITION;
if length(PARTITION) ~= size(X,1)
    disp('Incorrect number of coordinate pairs.');
    return
end
% h0 = image_figure;
hold on;
for i=1:length(PARTITION)
    if X(i,1)>=0
        plot(X(i,1),X(i,2),'Color',[.8 .8 .8]);
        text(X(i,1),X(i,2),num2str(PARTITION(i)));
    end
end
hold off;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function voronoiTessellation
% Tekee tulostiedostosta voronoi tessellaation.
h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
% waitALittle;
% [filename, pathname] = uigetfile('*.mat', 'Load mixture clustering');
% %load([pathname filename],'c');
% struct_array = load([pathname filename]);
% if isfield(struct_array,'c')  %Matlab versio
%     c = struct_array.c;
%     if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'PARTITION')  %Mideva versio
%     c = struct_array;
%     if ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% else
%     disp('Incorrect file format');
%     return;
% end
% 
% if ~isfield(c, 'pointers')
%     disp('Coordinate data missing from the result file');
%     return;
% end

pointers = c.pointers; vorPoints = c.vorPoints; vorCells = c.vorCells;
coordinates = c.coordinates; names = c.names;

if isequal(c.mixtureType, 'pop') || isequal(c.mixtureType, 'spatialPop')
    PARTITION = c.groupPartition;
else
    PARTITION = c.PARTITION;
end

talle = questdlg(['Do you want names to be visible in the colored ' ...
    'Voronoi tessellation?'], 'Names visible?', 'Yes', 'No', 'Yes');

if isequal(talle,'No')
    names = -1;
end
vorPlot(vorPoints, vorCells, PARTITION, pointers, coordinates, names);

%--------------------------------------------------------------------------

function localUncertainty
% Tekee tulostiedostosta kolmiulotteisen lokaalia epävarmuutta kuvaavan
% kuvan.

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
% waitALittle;
% [filename, pathname] = uigetfile('*.mat', 'Load mixture clustering');
% %load([pathname filename],'c');
% struct_array = load([pathname filename]);
% if isfield(struct_array,'c')  %Matlab versio
%     c = struct_array.c;
%     if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'PARTITION')  %Mideva versio
%     c = struct_array;
%     if ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% else
%     disp('Incorrect file format');
%     return;
% end
% 
% if ~isfield(c, 'pointers')
%     disp('Coordinate data missing from the result file');
%     return;
% end

pointers = c.pointers; vorPoints = c.vorPoints; vorCells = c.vorCells;
coordinates = c.coordinates; names = c.names;
varmuus = c.varmuus;

if isequal(c.mixtureType, 'pop') || isequal(c.mixtureType, 'spatialPop')
    PARTITION = c.groupPartition;
else
    PARTITION = c.PARTITION;
end

talle = questdlg('Do you want names to be visible in the plot?', ...
    'Names visible?', 'Yes', 'No', 'Yes');

if isequal(talle,'No')
    names = -1;
end

plotVarmuus(vorPoints, vorCells, pointers, varmuus, coordinates, ...
    PARTITION, names);



%--------------------------------------------------------------------------

function viewPartitio

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');

% waitALittle;
% [filename, pathname] = uigetfile('*.mat', 'Load mixture clustering');
% %load([pathname filename],'c');
% if (sum(filename)==0) || (sum(pathname)==0)
%     return;
% end
% struct_array = load([pathname filename]);
% if isfield(struct_array,'c')  %Matlab versio
%     c = struct_array.c;
%     if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'PARTITION')  %Mideva versio
%     c = struct_array;
%     if ~isfield(c,'rowsFromInd')
%         disp('Incorrect file format');
%         return
%     end
% else
%     disp('Incorrect file format');
%     return;
% end

if isequal(c.mixtureType, 'pop') || isequal(c.mixtureType, 'spatialPop')
    viewPopMixPartition(c.groupPartition, c.rows, c.popnames);
elseif isequal(c.mixtureType, 'trained')
    viewMixPartition(c.PARTITION, []);
else
    viewMixPartition(c.PARTITION, c.popnames);
end

function openHelpDoc
% s = fileparts(which('BAPS4manual.doc'));
% helpwin(s);
if strcmp(computer,'PCWIN')
    % s = fileparts(which('baps4.exe'));
    % winopen([s '\BAPS4manual.doc']);
    winopen('BAPS5manual.doc');
end

function openHelpHtml
% web http://www.rni.helsinki.fi/~jic/bapspage.html
% web('http://www.rni.helsinki.fi/~jic/bapspage.html','-browser')
% web http://www.rni.helsinki.fi/~jic/bapspage.html -new;
if strcmp(computer,'PCWIN')
    dos('start http://www.abo.fi/fak/mnf/mate/jc/software/baps.html'); % For the compiled version
end

function openAboutWindow
info{1}='';
info{2}='Bayesian Analysis of Population Structure (BAPS)';
info{3}='';
info{4}='Version 6.0';
info{5}='';
info{6}='Author: Jukka Corander, Pekka Marttinen, Jukka Siren, Jing Tang and Lu Cheng';
info{7}='';
info{8}='Copyright 2005-2012.  All Rights Reserved';
info{9}='';
info{10}='Please view the reference page when using as part of research';
info{11}='at http://www.helsinki.fi/bsg/software/BAPS';
info{12} ='';
helpdlg(info,'About');

%--------------------------------------------------------------------------

function viewAdmixture

% waitALittle;
% [filename, pathname] = uigetfile('*.mat', 'Load admixture results.');
% if (sum(filename)==0) || (sum(pathname)==0)
%     return;
% end
% %load([pathname filename],'c');
% struct_array = load([pathname filename]);
disp('---------------------------------------------------');
disp('Viewing the admixture result...');
h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
h0 = findobj('Tag','filename1_text');
filename = get(h0,'String');

% if isfield(struct_array,'c')  %Matlab versio
%     c = struct_array.c;
%     if ~isfield(c,'proportionsIt')
%         disp('*** ERROR: Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'proportionsIt')  %Mideva versio
%     c = struct_array;
%     if ~isfield(c,'proportionsIt')
%         disp('*** ERROR: Incorrect file format');
%         return
%     end
% else
%     disp('*** ERROR: Incorrect file format');
%     return;
% end

% mixtureType = c.mixtureType;
proportionsIt = c.proportionsIt; 
popnames = c.popnames; partition = c.PARTITION;
mixtureType = c.mixtureType;
% if strcmp(mixtureType,'linkage_mix') % For bacterial clustering data
% if isempty(popnames) || size(popnames,1)==size(partition,1)
%if  strcmp(mixtureType, 'admix')
    if isempty(popnames)
        ninds = size(partition,1);
        popnames=cell(ninds,2);
        for ind=1:ninds
            popnames{ind,1}=cellstr(num2str(ind));
        end        
        popnames(:,2)=num2cell((1:ninds)');
    end
 
    npops = c.npops; % all the clusters including outliers
    admixnpops = c.admixnpops;
    
    if ~isfield(c,'pvalue') % compatiable with old data
        disp('*** WARNING: pvalue is not found in the admixture result.');
        disp('*** WARNING: Old admixture file.');
        pvalue = ones(size(partition,1),1);
    else
        pvalue = c.pvalue;
    end
    view_admixture(proportionsIt,npops,admixnpops, ...
                         popnames,partition,pvalue,filename);
%else 
%    disp('*** ERROR: incorrect admixture data.');
    % put which variable as the input?
    % admixnpops = c.admixnpops;
%     npops = c.npops;
%     talle = questdlg(['Do you want individual names to be visible in the admixture ' ...
%         'result graphics?'], 'Names visible?', 'Yes', 'No', 'Yes');
%     if isequal(talle,'No')
%         viewPartition2(proportionsIt, [], npops, partition, filename);
%     else
%         viewPartition2(proportionsIt, popnames, npops, partition, filename);
%     end
% end

%--------------------------------------------------------------------------
function viewLoghood
    view_loglikelihood;
        
function viewEnergy
    view_energy;
            %--------------------------------------------------------------------------
 function viewGeneflow
    view_geneflow;
                
%--------------------------------------------------------------------------
function changeColor()
   h0 = findobj('Tag','base_figure');
   c = uisetcolor(h0,'Change color');
   h1 = findobj('Tag','datafile_text');
   h2 = findobj('Tag','outputfile_text');
   h3 = findobj('Tag','filename1_text');
   h4 = findobj('Tag','filename2_text');
   set(h1,'BackGroundColor',c);
   set(h2,'BackGroundColor',c);
   set(h3,'BackGroundColor',c);
   set(h4,'BackGroundColor',c);
   drawnow;
%-----------------------------------------------------------------------
function showmethebaps()
h0 = findobj('Tag','base_figure');
%setWindowOnTop(h0,'true')
goToDefault
h0 = findobj('Tag','load_menu');
set(h0,'UserData',[]);
%setWindowOnTop(h0,'false')


%-----------------------------------------------------------------------
function loadMixture

waitALittle;
[filename, pathname] = uigetfile('*.mat', 'Load mixture result');
%load([pathname filename],'c');
if (sum(filename)==0) || (sum(pathname)==0)
    return;
end
disp('---------------------------------------------------');
disp('In loading the mixture result...');
struct_array = load([pathname filename]);
if isfield(struct_array,'c')  %Matlab versio
    c = struct_array.c;
    if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd')
        disp('*** ERROR: Incorrect file format');
        return
    end
elseif isfield(struct_array,'PARTITION')  %Mideva versio
    c = struct_array;
    if ~isfield(c,'rowsFromInd')
        disp('*** ERROR: Incorrect file format');
        return
    end
else
    disp('*** ERROR: Incorrect file format');
    return;
end
% Save gathered information to 'mixture_menu's UserData:
h0 = findobj('Tag','load_menu');
set(h0,'UserData',c);
clear c;
%Set the name of the datafile to screen.
h1 = findobj('Tag','filename1_text');
if exist('pathname')
    filename = [pathname filename];
end
set(h1,'String',filename);
h1 = findobj('Tag','datafile_text');
set(h1,'String','Mixture result:');
disp('Mixture result loaded.');
goToMixtureAnalysis

%--------------------------------------------------------------------------
function goToMixtureAnalysis
set(findobj('Tag','graph_menu'), 'Enable','on');
set(findobj('Tag','partitio_menu'), 'Enable','on');
set(findobj('Tag','likelihood_menu'), 'Enable','on');
set(findobj('Tag','energy_menu'), 'Enable','on');
set(findobj('Tag','distances_menu'), 'Enable','on');
set(findobj('Tag','kl_menu'), 'Enable','on');
set(findobj('Tag','nei_menu'), 'Enable','on');
set(findobj('Tag','close_menu'), 'Enable','on');
set(findobj('Tag','phylogeny_menu'), 'Enable','on');
set(findobj('Tag','upgma_menu'), 'Enable','on');
set(findobj('Tag','nj_menu'), 'Enable','on');
set(findobj('Tag','geneflow_menu'), 'Enable','off');
set(findobj('Tag','admix_menu'), 'Enable','off');
set(findobj('Tag','mutationplot_menu'), 'Enable','on');



%--------------------------------------------------------------------------
function goToDefault
set(findobj('Tag','graph_menu'), 'Enable','off');
set(findobj('Tag','distances_menu'), 'Enable','off');




%--------------------------------------------------------------------------
function loadSpatial
% Tekee tulostiedostosta voronoi tessellaation.

waitALittle;
[filename, pathname] = uigetfile('*.mat', 'Load spatial mixture/admixture clustering');
%load([pathname filename],'c');
if (sum(filename)==0) || (sum(pathname)==0)
    return;
end
struct_array = load([pathname filename]);
disp('---------------------------------------------------');
disp('In loading the spatial mixture/admixture result...');
if isfield(struct_array,'c')  %Matlab versio
    c = struct_array.c;
    if ~isfield(c,'PARTITION') || ~isfield(c,'mixtureType')
        disp('*** ERROR: Incorrect file format');
        return
    end
    if ~strcmp(c.mixtureType,'spatial') && ~strcmp(c.mixtureType,'spatialPop')
        disp('*** ERROR: Incorrect file format');
        return
    end
elseif isfield(struct_array,'PARTITION')  %Mideva versio
    c = struct_array;
    if ~isfield(c,'rowsFromInd')
        disp('*** ERROR: Incorrect file format');
        return
    end
else
    disp('*** ERROR: Incorrect file format');
    return;
end

if ~isfield(c, 'pointers')
    disp('*** ERROR: Coordinate data missing from the result file');
    return;
end

% Save gathered information to 'spatialmixture_menu's UserData:
h0 = findobj('Tag','load_menu');
set(h0,'UserData',c);

%Set the name of the datafile to screen.
h1 = findobj('Tag','filename1_text');
if exist('pathname')
    filename = [pathname filename];
end
set(h1,'String',filename);
h1 = findobj('Tag','datafile_text');
if isfield(c,'admixnpops')
    set(h1,'String','Spatial AdMixture Result:');
    disp('Spatial admixture result loaded.');
    goToSpatialAdMixtureAnalysis
else
    set(h1,'String','Spatial Mixture Result:');
    disp('Spatial mixture result loaded.');
    goToSpatialMixtureAnalysis
end
clear c;


%--------------------------------------------------------------------------
function goToSpatialMixtureAnalysis
set(findobj('Tag','graph_menu'), 'Enable','on');
set(findobj('Tag','plot_coordinates_menu'), 'Enable','on');
set(findobj('Tag','voronoi_menu'), 'Enable','on');
set(findobj('Tag','varmuus_menu'), 'Enable','on');
set(findobj('Tag','distances_menu'), 'Enable','on');
set(findobj('Tag','kl_menu'), 'Enable','on');
set(findobj('Tag','nei_menu'), 'Enable','on');
set(findobj('Tag','close_menu'), 'Enable','on');
set(findobj('Tag','geneflow_menu'), 'Enable','off');
set(findobj('Tag','admix_menu'), 'Enable','off');
set(findobj('Tag','likelihood_menu'), 'Enable','on');
set(findobj('Tag','partitio_menu'), 'Enable','off');

%--------------------------------------------------------------------------
function goToSpatialAdMixtureAnalysis
set(findobj('Tag','graph_menu'), 'Enable','on');
set(findobj('Tag','plot_coordinates_menu'), 'Enable','on');
set(findobj('Tag','voronoi_menu'), 'Enable','on');
set(findobj('Tag','varmuus_menu'), 'Enable','on');
set(findobj('Tag','distances_menu'), 'Enable','on');
set(findobj('Tag','kl_menu'), 'Enable','on');
set(findobj('Tag','nei_menu'), 'Enable','on');
set(findobj('Tag','close_menu'), 'Enable','on');
set(findobj('Tag','geneflow_menu'), 'Enable','off');
set(findobj('Tag','admix_menu'), 'Enable','on');
set(findobj('Tag','likelihood_menu'), 'Enable','on');
set(findobj('Tag','partitio_menu'), 'Enable','off');
%--------------------------------------------------------------------------
function closeFile
h0 = findobj('Tag','load_menu');
set(h0,'UserData',[]);
h0 = findobj('Tag','datafile_text');
set(h0,'String','Data File:');
h0 = findobj('Tag','filename1_text');
set(h0,'String','');

set(findobj('Tag','close_menu'), 'Enable','off');
set(findobj('Tag','graph_menu'), 'Enable','off');
set(findobj('Tag','plot_coordinates_menu'), 'Enable','off');
set(findobj('Tag','partitio_menu'), 'Enable','off');
set(findobj('Tag','likelihood_menu'), 'Enable','off');
set(findobj('Tag','admix_menu'), 'Enable','off');
set(findobj('Tag','geneflow_menu'), 'Enable','off');
set(findobj('Tag','voronoi_menu'), 'Enable','off');
set(findobj('Tag','varmuus_menu'), 'Enable','off');

set(findobj('Tag','distances_menu'), 'Enable','off');
set(findobj('Tag','kl_menu'), 'Enable','off');
set(findobj('Tag','nei_menu'), 'Enable','off');


%-----------------------------------------------------------------------
function loadAdmixture
waitALittle;
[filename, pathname] = uigetfile('*.mat', 'Load admixture results.');
if (sum(filename)==0) || (sum(pathname)==0)
    return;
end
%load([pathname filename],'c');
disp('---------------------------------------------------');
disp('In loading the admixture result...');
struct_array = load([pathname filename]);

if isfield(struct_array,'c')  %Matlab versio
    c = struct_array.c;
    if ~isfield(c,'proportionsIt')
        disp('*** ERROR: Incorrect file format');
        return
    end
elseif isfield(struct_array,'proportionsIt')  %Mideva versio
    c = struct_array;
    if ~isfield(c,'proportionsIt')
        disp('*** ERROR: Incorrect file format');
        return
    end
elseif isfield(struct_array,'tietue')
    c = struct_array.tietue;
    if ~isfield(c,'proportionsIt')
        disp('*** ERROR: Incorrect file format');
        return
    end
else
    disp('*** ERROR: Incorrect file format');
    return;
end

% Save gathered information to 'mixture_menu's UserData:
h0 = findobj('Tag','load_menu');
set(h0,'UserData',c);
clear c;
%Set the name of the datafile to screen.
h1 = findobj('Tag','filename1_text');
if exist('pathname')
    filename = [pathname filename];
end
set(h1,'String',filename);
h1 = findobj('Tag','datafile_text');
set(h1,'String','Admixture result:');
disp('Admixture result loaded.');
goToAdmixtureAnalysis

%--------------------------------------------------------------------------
function goToAdmixtureAnalysis
set(findobj('Tag','graph_menu'), 'Enable','on');
set(findobj('Tag','admix_menu'), 'Enable','on');
set(findobj('Tag','geneflow_menu'), 'Enable','on');
% set(findobj('Tag','distances_menu'), 'Enable','on');
% set(findobj('Tag','kl_menu'), 'Enable','on');
% set(findobj('Tag','nei_menu'), 'Enable','on');
set(findobj('Tag','close_menu'), 'Enable','on');
% set(findobj('Tag','likelihood_menu'), 'Enable','on');


%--------------------------------------------------------------------------
function calculateDis(type)
if exist('baps4_output.baps','file')
    delete('baps4_output.baps')
else
    diary('baps4_output.baps')
end

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
npops = c.npops;
COUNTS = c.COUNTS;
adjprior = c.adjprior;
data = c.data;
partition = c.PARTITION;
clear c;

if npops > 1
    dist_mat = zeros(npops, npops);
    maxnoalle = size(COUNTS,1);
    nloci = size(COUNTS,2);
    d = zeros(maxnoalle, nloci, npops);
   
    
    switch type
        case 'KL'
            prior = adjprior;
            prior(find(prior==1))=0;
            nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
            prior(1,nollia)=1;
            for pop1 = 1:npops
                d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
                %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
            end

            ekarivi = num2str(npops);
            disp('--------------------------------------');
            disp('KL-divergence matrix in PHYLIP format:');
            disp('--------------------------------------');
            disp(ekarivi);
            for pop1 = 1:npops
                % rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
                for pop2 = 1:pop1-1
                    dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
                    div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
                    div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
                    div = (div12+div21)/2;
                    dist_mat(pop1,pop2) = div;
                end
            end

       case 'Nei'
           for pop1 = 1:npops
               d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))),maxnoalle,1);
               %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
           end

           ekarivi = num2str(npops);
           disp('--------------------------------------');
           disp('Nei-divergence matrix in PHYLIP format:');
           disp('--------------------------------------');
           disp(ekarivi);
           for pop1 = 1:npops
                % rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
                for pop2 = 1:pop1-1
                    dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
                    div1 = sum(sum(dist1.*dist2));
                    div2 = sqrt(sum(sum(dist1.^2)))*sqrt(sum(sum(dist2.^2)));
                    div = -log(div1/div2);
                    dist_mat(pop1,pop2) = div;
                end
           end    
        case 'Hamming'
           ekarivi = num2str(npops);
           disp('--------------------------------------');
           disp('Hamming distance matrix in PHYLIP format:');
           disp('--------------------------------------');
           disp(ekarivi);
            for pop1 = 1:npops
                for pop2 = 1:pop1-1
                    dist_mat(pop1,pop2) = hamming_dist(data(logical(partition==pop1),[1:end-1]),...
                                                  data(logical(partition==pop2),[1:end-1]));
                end
            end
    end
    
end

dist_mat = dist_mat + dist_mat'; % make it symmetric
for pop1 = 1:npops
    rivi = ['Cluster_' num2str(pop1) ' '];
    for pop2 = 1:npops
        rivi = [rivi kldiv2str(dist_mat(pop1,pop2)) ' '];
    end
    disp(rivi);
end
diary off

% ---------------------------------------------------------------------
% Save the result.
% Jing - 26.12.2005
talle = questdlg(['Do you want to save the distance matrix in PHYLIP format? '], ...
    'Save distance matrix?','Yes','No','Yes');
if isequal(talle,'Yes')
    %%%waitALittle;
    [filename, pathname] = uiputfile('*.txt','Save results as');

    if (sum(filename)==0) || (sum(pathname)==0)
        % Cancel was pressed
        return
    else % copy 'baps4_output.baps' into the text file with the same name.
        if exist('baps4_output.baps','file')
            copyfile('baps4_output.baps',[pathname filename])
            delete('baps4_output.baps')
        end
    end
else
    delete('baps4_output.baps')
end

% -------------------------------------------------------------------------
function mjono = kldiv2str(div)
mjono = '      ';
if abs(div)<100
    %Ei tarvita e-muotoa
    mjono(6) = num2str(rem(floor(div*1000),10));
    mjono(5) = num2str(rem(floor(div*100),10));
    mjono(4) = num2str(rem(floor(div*10),10));
    mjono(3) = '.';
    mjono(2) = num2str(rem(floor(div),10));
    arvo = rem(floor(div/10),10);
    if arvo>0
        mjono(1) = num2str(arvo);
    end
    
else
    suurinYks = floor(log10(div));
    mjono(6) = num2str(suurinYks);
    mjono(5) = 'e';
    mjono(4) = palautaYks(abs(div),suurinYks-1);
    mjono(3) = '.';
    mjono(2) = palautaYks(abs(div),suurinYks);
end

% -------------------------------------------------------------------------
function dist = hamming_dist(data1,data2)
[length1,nloci] = size(data1);
length2 = size(data2,1);
dist1 = 0;
for i = 1:length1
    dist2 = 0;
    for j = 1:length2
        dist2 = dist2 + sum(data1(i,:)~=data2(j,:))/nloci;
    end
    dist1 = dist1 + dist2/length2;
end
dist = dist1/length1;
%--------------------------------------------------------------------------
function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in?
% yks täytyy olla kokonaisluku, joka on 
% vähintään -1:n suuruinen. Pienemmill?
% luvuilla tapahtuu jokin pyöristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);




