function spatialPopMixture()
%Vaihtuvalla populaatioiden m‰‰r‰ll? priori 3:lla
%Samassa pisteess?olevien yksilˆiden t‰ytyy olla samasta populaatiosta.
%Toiminta pitk‰lti samanlainen kuin greedyPopMixiss?

base = findobj('Tag','base_figure'); % added by Lu Cheng, 11.11.2012

% check whether fixed k mode is selected
h0 = findobj('Tag','fixk_menu');
fixedK = get(h0, 'userdata');

if fixedK
    if ~(fixKWarning == 1) % call function fixKWarning
        return
    end
end

% output file name
OUTPUT_FILE = 'baps5_output.baps';  % also remember to update the file name in function WriteMixtureInfo

% check whether partition compare mode is selected
h1 = findobj('Tag','partitioncompare_menu');
partitionCompare = get(h1, 'userdata');

formatList = {'BAPS-format','FASTA-format', 'GenePop-format', 'Preprocessed data'};
formatChoice = menu('Specify the format of your data: ','BAPS-format','FASTA-format', 'GenePop-format', 'Preprocessed data');
if formatChoice==0
    return;
else
    input_type = formatList{formatChoice};
end

% input_type = questdlg('Specify the format of your data: ',...
%     'Specify Data Format', ...
%     'BAPS-format', 'GenePop-format', 'Preprocessed data', ...
%     'BAPS-format');

switch input_type
   
case 'BAPS-format'
    waitALittle;
    [filename1, pathname1] = uigetfile('*.txt', 'Load data in BAPS-format');
    if filename1==0
        return;
    end
    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname1 filename1]);
    end    
	
    data = load([pathname1 filename1]);
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    if (ninds==0) 
        disp('Incorrect Data-file.');
        return;    
    end
     
    waitALittle;
    [filename2,pathname2]=uigetfile('*.txt', 'Load group coordinates');
    if filename2==0
        return
    end
    
    coordinates = load([pathname2 filename2]);
    %viallinen = testaaKoordinaatit(ninds, coordinates);
    [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012
    if viallinen
        disp('Incorrect coordinates');
        return
    end
    
    inp = [filename1 ' & ' filename2];
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',inp);
    clear h0; clear inp;
    clear filename1; clear filename2; clear pathname1; clear pathname2;
        
    load_names = questdlg('Do you wish to specify the names of the groups?',...
        'Input group names?','Yes','No','Yes');
    if isequal(load_names,'Yes') 
        waitALittle;
        [filename, pathname] = uigetfile('*.txt', 'Load group names');
        popnames = initPopNames([pathname filename]);
        if (size(popnames,1)~=ninds)
            disp('Incorrect name-file.');
            popnames = [];
        end
    else
        popnames = [];
    end
    
    disp('Pre-processing the data. This may take several minutes.');
	   
    [cliques, separators, vorPoints, vorCells, pointers] = ...
        handleCoords(coordinates); 
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
    [Z,dist] = newGetDistances(data,rows);
    
    rowsFromInd = 0;  % Ei tiedet?
    	
    save_preproc = questdlg('Do you wish to save pre-processed data?',...
       'Save pre-processed data?',...
       'Yes','No','Yes');
    if isequal(save_preproc,'Yes');
        waitALittle;
        [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
        kokonimi = [pathname filename];
        c.data = data; c.rows = rows; c.alleleCodes = alleleCodes;
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.dist = dist; c.popnames = popnames; c.Z = Z;
        c.cliques = cliques; c.separators = separators;
        c.vorPoints = vorPoints; c.rowsFromInd = rowsFromInd;
        c.vorCells = vorCells; c.pointers = pointers;
        c.coordinates = coordinates;
%         save(kokonimi,'c');
        save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        clear c;
    end;
    
%%%%%%%%%%%%% added by Lu Cheng 11.11.2012 START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
case 'FASTA-format'

%     setWindowOnTop(base,'false')
    [filename1, pathname1] = uigetfile({'*.fasta';'*.*'}, 'Load data in FASTA-format');
    if filename1==0
        return;
    end
	
    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
    end
    %data = load([pathname1 filename1]);
    [~, seqs] = fastaread([pathname1 filename1]);
    seqs = seqs(:);
    alnMat = cell2mat(seqs);
    nSeq = length(seqs);
    clear seqs;
    
    cc = preprocAln(alnMat);
    
    setWindowOnTop(base,'false')
    [filename1_1,pathname1_1]=uigetfile('*.txt', 'Load group partition');
    if filename1_1==0
        return
    end
    groupPartition = load([pathname1_1 filename1_1]);
    if nSeq~=length(groupPartition)
        warning('Number of individuals inconsistent in %s (%d) and %s (%d).',...
            filename1,nSeq,filename1_1,length(groupPartition));
        return;
    else
        nPregroup = length(unique(groupPartition));
        ninds = nPregroup;
        assert(nPregroup==max(groupPartition));
        
        cc.nPregroup = nPregroup;
        cc.groupPartition = groupPartition;
        
        disp('Calculating distance matrix. This may take several minutes.');
        
        pgdist = nchoosek(nPregroup,2);
        tmpIndK=1;
        for i=1:nPregroup
            for j=i+1:nPregroup
                tmpIndsI = find(groupPartition==i);
                tmpIndsJ = find(groupPartition==j);
                tmpNI = length(tmpIndsI);
                tmpNJ = length(tmpIndsJ);
                tmpSum=0;
                for k=tmpIndsJ(:)'
                    tmp = alnMat(tmpIndsI,:)~=repmat(alnMat(k,:),tmpNI,1);
                    tmpSum = tmpSum+sum(tmp(:));
                end
                                
                pgdist(tmpIndK)=tmpSum/tmpNI/tmpNJ;
                tmpIndK = tmpIndK+1;
            end
        end
        clear tmp* alnMat
    end

    setWindowOnTop(base,'false')
    [filename2,pathname2]=uigetfile('*.txt', 'Load group coordinates');
    if filename2==0
        return
    end
    
    coordinates = load([pathname2 filename2]);
    %viallinen = testaaKoordinaatit(ninds, coordinates);
    [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012
    if viallinen
        disp('Incorrect coordinates');
        return
    end
    
    inp = [filename1 ' & ' filename1_1 ' & ' filename2];
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',inp);
    clear h0; clear inp;
    clear filename1; clear filename2; clear pathname1; clear pathname2;
    clear filename1_1 pathname1_1

    load_names = questdlg('Do you wish to specify the names of the groups?',...
        'Input group names?','Yes','No','Yes');
%     load_names = 'No';
    if isequal(load_names,'Yes') 
        waitALittle;
        [filename, pathname] = uigetfile('*.txt', 'Load group names');
        popnames = initPopNames([pathname filename]);
        if (size(popnames,1)~=ninds)
            disp('Incorrect name-file.');
            popnames = [];
        else
            popnames = popnames(:,1);
            
        end
    else
        popnames = [];
    end
    
    disp('Pre-processing the data. This may take several minutes.');
	   
    [cliques, separators, vorPoints, vorCells, pointers] = ...
        handleCoords(coordinates); 
    
    cc.locCliques = cliques;
    cc.locSeparators = separators;
    cc.popnames = popnames;
    cc.vorPoints = vorPoints; 
    cc.vorCells = vorCells;
    cc.pointers = pointers; 
    cc.coordinates = coordinates;
	format_type = 'FASTA';
    	
    save_preproc = questdlg('Do you wish to save pre-processed data?',...
       'Save pre-processed data?',...
       'Yes','No','Yes');
    if isequal(save_preproc,'Yes');
        waitALittle;
        [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
        kokonimi = [pathname filename];
        save(kokonimi,'cc','pgdist','format_type','groupPartition','-v7.3');
    end
    
    handlePopFastaCase(cc,groupPartition,pgdist);
    
    return;
    
%%%%%%%%%%%%%add by Lu Cheng 11.11.2012  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   
case 'GenePop-format'
    waitALittle;
    [filename1, pathname1] = uigetfile('*.txt', 'Load data in GenePop-format');
    if filename1==0
        return;
    end

    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname1 filename1]);
    end
    
    kunnossa = testaaGenePopData([pathname1 filename1]);
    if kunnossa==0
        return
    end
    [data,popnames]=lueGenePopData([pathname1 filename1]);
    
    waitALittle;	   
    [filename2,pathname2]=uigetfile('*.txt', 'Load group coordinates');
    if filename2==0
        return
    end
    
    ninds = max(data(:,end));    
    coordinates = load([pathname2 filename2]);
    %viallinen = testaaKoordinaatit(ninds, coordinates);
    [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012
    
    if viallinen
        disp('Incorrect coordinates');
        return
    end
    
    inp = [filename1 ' & ' filename2];
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',inp);
    clear h0; clear inp;
    clear filename1; clear filename2; clear pathname1; clear pathname2;
    
    disp('Pre-processing the data. This may take several minutes.');
    
    [cliques, separators, vorPoints, vorCells, pointers] = ...
        handleCoords(coordinates);    
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
    [Z,dist] = newGetDistances(data,rows);
    
    rowsFromInd = 2; %Tiedet‰‰n
    
    save_preproc = questdlg('Do you wish to save pre-processed data?',...
       'Save pre-processed data?',...
       'Yes','No','Yes');
    if isequal(save_preproc,'Yes');
        waitALittle;
        [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
        kokonimi = [pathname filename];
        c.data = data; c.rows = rows; c.alleleCodes = alleleCodes;
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.dist = dist; c.popnames = popnames; c.Z = Z;
        c.cliques = cliques; c.separators = separators;       
        c.vorPoints = vorPoints; c.rowsFromInd = rowsFromInd;
        c.vorCells = vorCells; c.pointers = pointers; 
        c.coordinates = coordinates;
%         save(kokonimi,'c');
        save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        clear c;
    end;      
   
case 'Preprocessed data'
    waitALittle;
    [filename, pathname] = uigetfile('*.mat', 'Load pre-processed data');
    if filename==0
        return;
    end
    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
    end    
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;

    struct_array = load([pathname filename]);
    
    if isfield(struct_array,'format_type') && strcmp(struct_array.format_type,'FASTA')
        handlePopFastaCase(struct_array.cc,struct_array.groupPartition,struct_array.pgdist);
        return;
    end
    
    if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        if ~isfield(c,'dist')
            disp('Incorrect file format');
            return
        end
    elseif isfield(struct_array,'dist')  %Mideva versio
        c = struct_array;
    else
        disp('Incorrect file format');
        return;
    end
    data = double(c.data); rows = c.rows; alleleCodes = c.alleleCodes;
    noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
    dist = c.dist; popnames = c.popnames; Z = c.Z; rowsFromInd = c.rowsFromInd;
    
    if isfield(c, 'cliques')
        cliques = c.cliques; separators = c.separators;
        vorPoints = c.vorPoints; vorCells = c.vorCells;
        pointers = c.pointers; coordinates = c.coordinates;
        clear c;
    else
        load_coord = questdlg(['The data file did not contain ',...
            'coordinate information. Do you wish to load coordinates?'], ...
            'Load coordinates?',...
            'Yes','No','Yes');
        if isequal(load_coord, 'No')
            return
        end
        waitALittle;
        [filename2,pathname2]=uigetfile('*.txt', 'Load group coordinates');
        if filename2==0
            return
        end

        ninds = max(data(:,end));
        coordinates = load([pathname2 filename2]);
        %viallinen = testaaKoordinaatit(ninds, coordinates);
        [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012

        if viallinen
            disp('Incorrect coordinates');
            return
        end
        inp = [filename ' & ' filename2];
        h0 = findobj('Tag','filename1_text');
        set(h0,'String',inp);
        clear h0; clear inp;
        clear filename; clear filename2; clear pathname; clear pathname2;
        
        disp('Pre-processing the data. This may take several minutes.');    
        [cliques, separators, vorPoints, vorCells, pointers] = ...
            handleCoords(coordinates);
        save_preproc = questdlg('Do you wish to save pre-processed data?',...
            'Save pre-processed data?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
            kokonimi = [pathname filename];
            c.cliques = cliques; c.separators = separators;
            c.vorPoints = vorPoints; c.vorCells = vorCells;
            c.pointers = pointers; c.coordinates = coordinates;
%             save(kokonimi,'c');
            save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
            clear c;
        end; 
    end
    otherwise
        return
end

global PARTITION; global COUNTS;
global SUMCOUNTS; %global POP_LOGML;
global SEPCOUNTS; global CLIQCOUNTS;
clearGlobalVars;

c.data=data; c.alleleCodes = alleleCodes;
c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
c.dist=dist; c.Z=Z; c.rowsFromInd = rowsFromInd;
c.cliques = cliques; c.separators = separators;
c.rows = rows;

% partition compare
if ~isempty(partitionCompare)
    nsamplingunits = size(rows,1);
    data = data(:,1:end-1);
    partitions = partitionCompare.partitions;
    npartitions = size(partitions,2);
    partitionLogml = zeros(1,npartitions);
    for i = 1:npartitions
        % number of unique partition lables
        npops = length(unique(partitions(:,i)));
        try
            partitionInd = zeros(rows(end),1);
            partitionSample = partitions(:,i);
            for j = 1:nsamplingunits
                partitionInd([c.rows(j,1):c.rows(j,2)]) = partitionSample(j);
            end
            [sumcounts, counts] = initialCounts(partitionInd, data, npops, noalle);
            [cliqcounts, sepcounts] = computeCounts(cliques, separators, npops, partitionSample);
            partitionLogml(i) = ...
                  computeLogml(adjprior, priorTerm, cliqcounts, sepcounts, counts, sumcounts);      
        catch
            disp('*** ERROR: unmatched data.');
            return
        end
    end
    % return the logml result
    partitionCompare.logmls = partitionLogml;
    set(h1, 'userdata', partitionCompare);
    return
end

if fixedK
    [logml, npops, partitionSummary]=spatialMix_fixK(c);
else    
    [logml, npops, partitionSummary]=spatialMix(c);
end

if logml==1
    return
end

lastCol = data(:,end);
data = data(:,1:end-1);

h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
h0 = findobj('Tag','filename2_text');
outp = get(h0,'String');
[varmuus,changesInLogml] = writeMixtureInfo(logml, rows, data, adjprior, priorTerm, ...
    outp,inp,partitionSummary, popnames, cliques, separators, fixedK);

%checkLogml(priorTerm, adjprior, cliques, separators);

viewPopMixPartition(PARTITION, rows, popnames);

if isequal(popnames, [])
    names = pointers;
else
   %Etsit‰‰n voronoi-soluja vastaavat nimet.
   names = cell(size(pointers));
   indices = 1:length(popnames);
   for i = 1:length(pointers)
       inds = pointers{i};
       namesInCell = [];
       for j = 1:length(inds)
           ind = inds(j);
           I = find(indices > ind);
           if isempty(I)
               nameIndex = indices(end);
           else
               nameIndex = min(I) -1;
           end
           name = popnames{nameIndex};
           namesInCell = [namesInCell name];
       end
       names{i} = namesInCell;
   end
end
vorPlot(vorPoints, vorCells, PARTITION, pointers, coordinates, names);

talle = questdlg(['Do you want to save the mixture populations ' ...
    'so that you can use them later in admixture analysis or plot ' ...
    'additional images?'], ...
    'Save results?','Yes','No','Yes');
if isequal(talle,'Yes')
    waitALittle;
    [filename, pathname] = uiputfile('*.mat','Save results as');
    
    if (filename == 0) & (pathname == 0)
        % Cancel was pressed
        return
    else % copy OUTPUT_FILE into the text file with the same name.
        if exist(OUTPUT_FILE,'file')
            copyfile(OUTPUT_FILE,[pathname filename '.txt'])
            delete(OUTPUT_FILE)
        end
    end
    
    if rowsFromInd==0 
        %K‰ytettiin BAPS-formaattia, eik?rowsFromInd ole tunnettu.
        [popnames, rowsFromInd] = findOutRowsFromInd(popnames, rows);
    end
    
    groupPartition = PARTITION;
    
    fiksaaPartitioYksiloTasolle(rows, rowsFromInd);

    c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
    c.alleleCodes = alleleCodes; c.adjprior = adjprior; c.popnames = popnames;
    c.rowsFromInd = rowsFromInd; c.data = data; c.npops = npops;
    c.noalle = noalle; c.groupPartition = groupPartition;
    c.pointers = pointers; c.vorPoints = vorPoints; c.vorCells = vorCells;
    c.coordinates = coordinates; c.names = names; c.varmuus = varmuus; 
    c.rows = rows; c.mixtureType = 'spatialPop'; 
    c.logml = logml; c.changesInLogml = changesInLogml;
    
    %  added by Lu Cheng, 05.12.2012    
    tmpFile = [pathname filename '.mapfile.txt'];
    fid = fopen(tmpFile,'w+');
    fprintf(fid,'GroupLabel\tLatitude\tLongitude\tDescription\tLabel\n');
    for i=1:max(lastCol)
        fprintf(fid,'%d\t%.10f\t%.10f\t%d_%d\t%d\n',i,coordinates(i,1),coordinates(i,2),...
                i,PARTITION(i),PARTITION(i));
    end
    fclose(fid);
    
%     save([pathname filename], 'c');
    save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
else
    if exist(OUTPUT_FILE,'file')
    delete(OUTPUT_FILE)
    end
end





%-------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------

function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
%global POP_LOGML; POP_LOGML = [];
global SEPCOUNTS; SEPCOUNTS = [];
global CLIQCOUNTS; CLIQCOUNTS = [];

%-------------------------------------------------------------------------------------


function rows = computeRows(rowsFromInd, inds, ninds)
% On annettu yksilˆt inds. Funktio palauttaa vektorin, joka
% sis‰lt‰‰ niiden rivien numerot, jotka sis‰lt‰v‰t yksilˆiden
% dataa.

rows = inds(:, ones(1,rowsFromInd));
rows = rows*rowsFromInd;
miinus = repmat(rowsFromInd-1 : -1 : 0, [ninds 1]);
rows = rows - miinus;
rows = reshape(rows', [1,rowsFromInd*ninds]);


%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedet‰‰n, ett?annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss?ei viel?ole
% annettua logml arvoa, niin lis‰t‰‰n worstIndex:in kohtaan uusi logml ja
% nykyist?partitiota vastaava nclusters:in arvo. Muutoin ei tehd?mit‰‰n.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lˆydetty partitio ei ole viel?kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end


%--------------------------------------------------------------------------


function [suurin, i2] = arvoSeuraavaTila(muutokset, logml)
% Suorittaa yksilˆn seuraavan tilan arvonnan

y = logml + muutokset;  % siirron j‰lkeiset logml:t
y = y - max(y);
y = exp(y);
summa = sum(y);
y = y/summa;
y = cumsum(y);

i2 = rand_disc(y);   % uusi kori
suurin = muutokset(i2);


%--------------------------------------------------------------------------------------


function svar=rand_disc(CDF)
%returns an index of a value from a discrete distribution using inversion method
slump=rand;
har=find(CDF>slump);
svar=har(1);


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, diffInCounts, ...
    cliques, separators, adjprior, priorTerm)
% Suorittaa globaalien muuttujien muutokset, kun yksil?ind
% siirret‰‰n koriin i2.

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
global CLIQCOUNTS;
global SEPCOUNTS;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

diffInCliqCounts = computeDiffInCliqCounts(cliques, ind);
diffInSepCounts = computeDiffInCliqCounts(separators, ind);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) - diffInCliqCounts;
CLIQCOUNTS(:,i2) = CLIQCOUNTS(:,i2) + diffInCliqCounts;
SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) - diffInSepCounts;
SEPCOUNTS(:,i2) = SEPCOUNTS(:,i2) + diffInSepCounts;

%POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%---------------------------------------------------------------------------------


function updateGlobalVariables2(i1, i2, diffInCounts, ...
    cliques, separators, adjprior, priorTerm);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilˆt siirret‰‰n koriin i2.

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
%global POP_LOGML;
global CLIQCOUNTS;
global SEPCOUNTS;

inds = find(PARTITION==i1);
PARTITION(inds) = i2;

diffInCliqCounts = CLIQCOUNTS(:,i1);
diffInSepCounts = SEPCOUNTS(:,i1);


COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

CLIQCOUNTS(:,i1) = 0;
CLIQCOUNTS(:,i2) = CLIQCOUNTS(:,i2) + diffInCliqCounts;
SEPCOUNTS(:,i1) = 0;
SEPCOUNTS(:,i2) = SEPCOUNTS(:,i2) + diffInSepCounts;


%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, ...
    adjprior, priorTerm, i2, cliques, separators);
% Suorittaa globaalien muuttujien p‰ivitykset, kun yksilˆt 'muuttuvat'
% siirret‰‰n koriin i2. Ennen siirtoa yksilˆiden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;      global CLIQCOUNTS;
global SUMCOUNTS;   global SEPCOUNTS;
%global POP_LOGML;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

diffInCliqCounts = computeDiffInCliqCounts(cliques, muuttuvat);
diffInSepCounts = computeDiffInCliqCounts(separators, muuttuvat);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) - diffInCliqCounts;
CLIQCOUNTS(:,i2) = CLIQCOUNTS(:,i2) + diffInCliqCounts;
SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) - diffInSepCounts;
SEPCOUNTS(:,i2) = SEPCOUNTS(:,i2) + diffInSepCounts;

%POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%----------------------------------------------------------------------


function inds = returnInOrder(inds, pop, globalRows, data, ...
    adjprior, priorTerm)
% Palauttaa yksilˆt j‰rjestyksess?siten, ett?ensimm‰isen?on
% se, jonka poistaminen populaatiosta pop nostaisi logml:n
% arvoa eniten.

global COUNTS;      global SUMCOUNTS;
ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind =inds(i);
    rows = globalRows(i,1):globalRows(i,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop) = COUNTS(:,:,pop)-diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)-diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
    COUNTS(:,:,pop) = COUNTS(:,:,pop)+diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)+diffInSumCounts;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset(ind, globalRows, ...
    data, adjprior, priorTerm, logml, cliques, separators)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li yksilˆt inds siirret‰‰n koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lis‰tt‰v?
% COUNTS:in siivuun i2, mik‰li muutos toteutetaan.
% Huom! Laskee muutoksen vain yhdelle tyhj‰lle populaatiolle, muiille
% tyhjille tulee muutokseksi 0.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   %global POP_LOGML;
global CLIQCOUNTS;  global SEPCOUNTS;

npops = size(COUNTS,3);
muutokset = zeros(npops,1);

counts = COUNTS;
sumcounts = SUMCOUNTS;

[emptyPop, pops] = findEmptyPop(npops);

i1 = PARTITION(ind);

i2 = [pops(find(pops~=i1))];
if emptyPop > 0
    i2 =[i2 emptyPop];
end

i2 = sort(i2);

rows = globalRows(ind,1):globalRows(ind,2);
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

diffInCliqCounts = computeDiffInCliqCounts(cliques, ind);
diffInSepCounts = computeDiffInCliqCounts(separators, ind);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) - diffInCliqCounts;
SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) - diffInSepCounts;

for i=i2
    CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) + diffInCliqCounts;
    SEPCOUNTS(:,i) = SEPCOUNTS(:,i) + diffInSepCounts;
    COUNTS(:,:,i) = COUNTS(:,:,i) + diffInCounts;
    SUMCOUNTS(i,:) = SUMCOUNTS(i,:) + diffInSumCounts;   
    
    muutokset(i) = computeLogml(adjprior, priorTerm) - logml;
    
    CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) - diffInCliqCounts;
    SEPCOUNTS(:,i) = SEPCOUNTS(:,i) - diffInSepCounts;    
    COUNTS(:,:,i) = COUNTS(:,:,i) - diffInCounts;
    SUMCOUNTS(i,:) = SUMCOUNTS(i,:) - diffInSumCounts;
end

COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;
CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) + diffInCliqCounts;
SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) + diffInSepCounts;

% Asetetaan muillekin tyhjille populaatioille sama muutos, kuin
% emptyPop:lle

if emptyPop > 0
    empties = mysetdiff((1:npops), [i2 i1]);
    muutokset(empties) = muutokset(emptyPop);
end

COUNTS = counts;
SUMCOUNTS = sumcounts;

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2(i1, globalRows, ...
    data, adjprior, priorTerm, logml, cliques, separators);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li korin i1 kaikki yksilˆt siirret‰‰n
% koriin i. 
% Laskee muutokset vain yhdelle tyhj‰lle populaatiolle, muille tulee
% muutokseksi 0.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global CLIQCOUNTS;  global SEPCOUNTS;

npops = size(COUNTS,3);
muutokset = zeros(npops,1);

[emptyPop, pops] = findEmptyPop(npops);

i2 = [pops(find(pops~=i1))];
if emptyPop > 0
    i2 =[i2 emptyPop];
end

inds = find(PARTITION == i1);
ninds = length(inds);

rows = [];
for i = 1:ninds
    rows = [rows globalRows(inds(i),1):globalRows(inds(i),2)];
end
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);
diffInCliqCounts = computeDiffInCliqCounts(cliques, inds);
diffInSepCounts = computeDiffInCliqCounts(separators, inds);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
CLIQCOUNTS(:,i1) = 0;
SEPCOUNTS(:,i1) = 0;

for i=i2
    CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) + diffInCliqCounts;
    SEPCOUNTS(:,i) = SEPCOUNTS(:,i) + diffInSepCounts;
    COUNTS(:,:,i) = COUNTS(:,:,i) + diffInCounts;
    SUMCOUNTS(i,:) = SUMCOUNTS(i,:) + diffInSumCounts;   
    
    muutokset(i) = computeLogml(adjprior, priorTerm) - logml;
    
    CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) - diffInCliqCounts;
    SEPCOUNTS(:,i) = SEPCOUNTS(:,i) - diffInSepCounts;    
    COUNTS(:,:,i) = COUNTS(:,:,i) - diffInCounts;
    SUMCOUNTS(i,:) = SUMCOUNTS(i,:) - diffInSumCounts;
end
    
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;
CLIQCOUNTS(:,i1) = diffInCliqCounts;
SEPCOUNTS(:,i1) = diffInSepCounts;



%------------------------------------------------------------------------------------

function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
    data, adjprior, priorTerm, i1, logml, cliques, separators)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mik?olisi muutos logml:ss? jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirret‰‰n koriin j.
% Laskee vain yhden tyhj‰n populaation, muita kohden muutokseksi j‰‰ 0.


global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global CLIQCOUNTS;  global SEPCOUNTS;

npops = size(COUNTS,3);
npops2 = length(unique(T2));
muutokset = zeros(npops2, npops);

for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
        rows = [];
        for i = 1:ninds
            ind = inds(i);
            rows = [rows; (globalRows(ind,1):globalRows(ind,2))'];
        end
        diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
        diffInSumCounts = sum(diffInCounts);
        diffInCliqCounts = computeDiffInCliqCounts(cliques, inds);
        diffInSepCounts = computeDiffInCliqCounts(separators, inds);
        
        COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
        CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) - diffInCliqCounts;
        SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) - diffInSepCounts;
        
        [emptyPop, pops] = findEmptyPop(npops);
        i2 = [pops(find(pops~=i1))];
        if emptyPop > 0
            i2 =[i2 emptyPop];
        end

        for i = i2
            CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) + diffInCliqCounts;
            SEPCOUNTS(:,i) = SEPCOUNTS(:,i) + diffInSepCounts;
            COUNTS(:,:,i) = COUNTS(:,:,i) + diffInCounts;
            SUMCOUNTS(i,:) = SUMCOUNTS(i,:) + diffInSumCounts;
                        
            muutokset(pop2,i) = computeLogml(adjprior, priorTerm) - logml;
                        
            CLIQCOUNTS(:,i) = CLIQCOUNTS(:,i) - diffInCliqCounts;
            SEPCOUNTS(:,i) = SEPCOUNTS(:,i) - diffInSepCounts;
            COUNTS(:,:,i) = COUNTS(:,:,i) - diffInCounts;
            SUMCOUNTS(i,:) = SUMCOUNTS(i,:) - diffInSumCounts;
        end

        COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;
        CLIQCOUNTS(:,i1) = CLIQCOUNTS(:,i1) + diffInCliqCounts;
        SEPCOUNTS(:,i1) = SEPCOUNTS(:,i1) + diffInSepCounts;
    end    
end

%--------------------------------------------------------------------------
function muutokset = laskeMuutokset5(inds, globalRows, data, ...
    adjprior, priorTerm, logml, cliques, separators, i1, i2)

% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li yksil?i vaihtaisi koria i1:n ja i2:n v‰lill?
    
global COUNTS;      global SUMCOUNTS;
global PARTITION;   
global CLIQCOUNTS;  global SEPCOUNTS;

ninds = length(inds);
muutokset = zeros(ninds,1);

for i = 1:ninds
    ind = inds(i);
    
    rows = globalRows(ind,1):globalRows(ind,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);
    
    if PARTITION(ind)==i1
        pop1 = i1;  %mist?
        pop2 = i2;  %mihin
    else
        pop1 = i2;
        pop2 = i1;
    end
        
    diffInCliqCounts = computeDiffInCliqCounts(cliques, ind);
    diffInSepCounts = computeDiffInCliqCounts(separators, ind);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)-diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)-diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)+diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)+diffInSumCounts;
    
    CLIQCOUNTS(:,pop1) = CLIQCOUNTS(:,pop1) - diffInCliqCounts;
    CLIQCOUNTS(:,pop2) = CLIQCOUNTS(:,pop2) + diffInCliqCounts;
    SEPCOUNTS(:,pop1) = SEPCOUNTS(:,pop1) - diffInSepCounts;
    SEPCOUNTS(:,pop2) = SEPCOUNTS(:,pop2) + diffInSepCounts;
    
    muutokset(i) = computeLogml(adjprior, priorTerm) - logml;
    
    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)+diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)+diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)-diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)-diffInSumCounts; 
    
    CLIQCOUNTS(:,pop1) = CLIQCOUNTS(:,pop1) + diffInCliqCounts;
    CLIQCOUNTS(:,pop2) = CLIQCOUNTS(:,pop2) - diffInCliqCounts;
    SEPCOUNTS(:,pop1) = SEPCOUNTS(:,pop1) + diffInSepCounts;
    SEPCOUNTS(:,pop2) = SEPCOUNTS(:,pop2) - diffInSepCounts;
end
    
%--------------------------------------------------------------------------

function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukum‰‰r‰t (vastaavasti kuin COUNTS:issa), jotka ovat data:n 
% riveill?rows.

diffInCounts = zeros(max_noalle, nloci);
for i=rows
    row = data(i,:);
    notEmpty = find(row>=0);
    
    if length(notEmpty)>0
        diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
            diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
    end
end    



%------------------------------------------------------------------------------------


function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on m‰‰ritelty pops-muuttujalla.

global COUNTS;
global SUMCOUNTS;
x = size(COUNTS,1);
y = size(COUNTS,2);
z = length(pops);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;

%------------------------------------------------------------------------------------

function npops = poistaTyhjatPopulaatiot(npops)
% Poistaa tyhjentyneet populaatiot COUNTS:ista ja 
% SUMCOUNTS:ista. P‰ivitt‰‰ npops:in ja PARTITION:in.

global COUNTS;
global SUMCOUNTS;
global PARTITION;
global CLIQCOUNTS;
global SEPCOUNTS;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);
CLIQCOUNTS = CLIQCOUNTS(:,notEmpty);
SEPCOUNTS = SEPCOUNTS(:,notEmpty);

for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
end
npops = length(notEmpty);


%----------------------------------------------------------------------------------
%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

function initial_partition=admixture_initialization(data_matrix,nclusters,Z)
size_data=size(data_matrix);
nloci=size_data(2)-1;
n=max(data_matrix(:,end));
T=cluster_own(Z,nclusters);
initial_partition=zeros(size_data(1),1);
for i=1:n
    kori=T(i);
    here=find(data_matrix(:,end)==i);
    for j=1:length(here)
        initial_partition(here(j),1)=kori;
    end
end

function T = cluster_own(Z,nclust)
true=logical(1);
false=logical(0);
maxclust = nclust;
% Start of algorithm
m = size(Z,1)+1;
T = zeros(m,1);
   % maximum number of clusters based on inconsistency
   if m <= maxclust
      T = (1:m)';
   elseif maxclust==1
      T = ones(m,1);
   else
      clsnum = 1;
      for k = (m-maxclust+1):(m-1)
         i = Z(k,1); % left tree
         if i <= m % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
         elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
         end
         i = Z(k,2); % right tree
         if i <= m  % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
         elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
         end
      end
   end
   
function T = clusternum(X, T, k, c)
m = size(X,1)+1;
while(~isempty(k))
   % Get the children of nodes at this level
   children = X(k,1:2);
   children = children(:);

   % Assign this node number to leaf children
   t = (children<=m);
   T(children(t)) = c;
   
   % Move to next level
   k = children(~t) - m;
end


%---------------------------------------------------------------------------------------


function [newData, rows, alleleCodes, noalle, adjprior, priorTerm] = ...
    handleData(raw_data)
% Alkuper‰isen datan viimeinen sarake kertoo, milt?yksilˆlt?
% kyseinen rivi on per‰isin. Funktio tutkii ensin, ett?montako
% rivi?maksimissaan on per‰isin yhdelt?yksilˆlt? jolloin saadaan
% tiet‰‰ onko kyseess?haploidi, diploidi jne... T‰m‰n j‰lkeen funktio
% lis‰‰ tyhji?rivej?niille yksilˆille, joilta on per‰isin v‰hemm‰n
% rivej?kuin maksimim‰‰r?
%   Mik‰li jonkin alleelin koodi on =0, funktio muuttaa t‰m‰n alleelin
% koodi pienimm‰ksi koodiksi, joka isompi kuin mik‰‰n k‰ytˆss?oleva koodi.
% T‰m‰n j‰lkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
% koodit saavat arvoja v‰lill?1,...,noalle(j).
%
% Muutettu vastaamaan greedyPopMixin handlePopDataa.

data = raw_data;
nloci=size(raw_data,2)-1;

dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
   isoinAlleeli = max(max(dataApu));
   dataApu(nollat) = isoinAlleeli+1;
   data(:,1:nloci) = dataApu;
end
dataApu = []; nollat = []; isoinAlleeli = [];

noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
end
alleleCodes = zeros(max(noalle),nloci);
for i=1:nloci
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
end

for loc = 1:nloci
    for all = 1:noalle(loc)
        data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;

nind = max(data(:,end));
%rows = cell(nind,1);
rows = zeros(nind,2);
for i=1:nind
    rivit = find(data(:,end)==i)';
    rows(i,1) = min(rivit);
    rows(i,2) = max(rivit);
end
newData = data;

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j)); 
end


%----------------------------------------------------------------------------------------

function [Z, dist] = newGetDistances(data, initRows)

ninds = size(initRows,1);
nloci = size(data,2)-1;
riviLkm = nchoosek(ninds,2);

empties = find(data<0);
data(empties)=0;
data = uint8(data);   % max(noalle) oltava <256

pariTaulu = zeros(riviLkm,2);
aPointer=1;
for a=1:ninds-1
    pariTaulu(aPointer:aPointer+ninds-1-a,1) = ones(ninds-a,1)*a;
    pariTaulu(aPointer:aPointer+ninds-1-a,2) = (a+1:ninds)';
    aPointer = aPointer+ninds-a;
end 

%eka = pariTaulu(:,ones(1,rowsFromInd));
%eka = eka * rowsFromInd;
%miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
%eka = eka - miinus;

koot = initRows(:,2) - initRows(:,1);
maxSize = max(koot) + 1;

rows = zeros(ninds, maxSize);

for i=1:ninds
    apu = initRows(i,1):initRows(i,2);
    rows(i, 1:length(apu)) = apu;
end
eka = zeros(riviLkm, maxSize);
toka = zeros(riviLkm, maxSize);

for i = 1:riviLkm
    eka(i, :) = rows(pariTaulu(i, 1), :);
    toka(i, :) = rows(pariTaulu(i,2), :);
end

%eka = uint16(eka);
%toka = uint16(toka);

summa = zeros(riviLkm,1);
vertailuja = zeros(riviLkm,1);

clear pariTaulu; clear miinus;

x = zeros(size(eka));    x = uint8(x);
y = zeros(size(toka));   y = uint8(y);

for j=1:nloci;
    
    for k=1:maxSize
        I = find(eka(:,k)>0);
        x(I,k) = data(eka(I,k),j);
        I = find(toka(:,k)>0);
        y(I,k) = data(toka(I,k),j);
    end
        
    for a=1:maxSize
        for b=1:maxSize
            vertailutNyt = double(x(:,a)>0 & y(:,b)>0);
            vertailuja = vertailuja + vertailutNyt;
            lisays = (x(:,a)~=y(:,b) & vertailutNyt);
            summa = summa+double(lisays);
        end
    end
end

clear x;    clear y;   clear vertailutNyt;
nollat = find(vertailuja==0);
dist = zeros(length(vertailuja),1);
dist(nollat) = 1;
muut = find(vertailuja>0);
dist(muut) = summa(muut)./vertailuja(muut);
clear summa; clear vertailuja;

Z = linkage(dist');

%----------------------------------------------------------------------------------------


function [Z, distances]=getDistances(data_matrix,nclusters)

%finds initial admixture clustering solution with nclusters clusters, uses simple mean Hamming distance
%gives partition in 8-bit format
%allocates all alleles of a single individual into the same basket
%data_matrix contains #Loci+1 columns, last column indicate whose alleles are placed in each row,
%i.e. ranges from 1 to #individuals. For diploids there are 2 rows per individual, for haploids only a single row
%missing values are indicated by zeros in the partition and by negative integers in the data_matrix.

size_data=size(data_matrix);
nloci=size_data(2)-1;
n=max(data_matrix(:,end));
distances=zeros(nchoosek(n,2),1);
pointer=1;
for i=1:n-1
    i_data=data_matrix(find(data_matrix(:,end)==i),1:nloci);
    for j=i+1:n
        d_ij=0;
        j_data=data_matrix(find(data_matrix(:,end)==j),1:nloci);
        vertailuja = 0;
        for k=1:size(i_data,1)
            for l=1:size(j_data,1)
                here_i=find(i_data(k,:)>=0);
                here_j=find(j_data(l,:)>=0);
                here_joint=intersect(here_i,here_j);
                vertailuja = vertailuja + length(here_joint);
                d_ij = d_ij + length(find(i_data(k,here_joint)~=j_data(l,here_joint)));
            end
        end
        d_ij = d_ij / vertailuja;
        distances(pointer)=d_ij;
        pointer=pointer+1;
    end
end

Z=linkage(distances');



%----------------------------------------------------------------------------------------


function Z = linkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 | m ~= fix(m)
  error('The first input has to match the output of the PDIST function in size.');   
end
if nargin == 1 % set default switch to be 'co' 
   method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
monotonic = 1;
Z = zeros(m-1,3); % allocate the output matrix.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n. 
R = 1:n;
for s = 1:(n-1)
   X = Y;
   [v, k] = min(X);
   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;
   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A   
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];
   
   switch method
   case 'si' %single linkage
      Y(I) = min(Y(I),Y(J));
   case 'av' % average linkage
      Y(I) = Y(I) + Y(J);
   case 'co' %complete linkage
      Y(I) = max(Y(I),Y(J));
   case 'ce' % centroid linkage
      K = N(R(i))+N(R(j));
      Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
   case 'wa'
      Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
	  N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
   end
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = []; % no need for the cluster information about j.
   
   % update m, N, R
   m = m-1; 
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n); 
end


%-----------------------------------------------------------------------------------


function popnames = initPopNames(nameFile)

fid = fopen(nameFile);
if fid == -1
    %File didn't exist
    msgbox('Loading of the population names was unsuccessful', ...
        'Error', 'error');
    return;
end;
line = fgetl(fid);
counter = 1;
while (line ~= -1) & ~isempty(line)
    names{counter} = line;
    line = fgetl(fid);
    counter = counter + 1;
end;
fclose(fid);

popnames = cell(length(names), 2);
for i = 1:length(names)
    popnames{i,1} = names(i);
    popnames{i,2} = 0;
end


%-----------------------------------------------------------------------------------
% Laskee arvot cliqcounts:lle ja sepcounts:lle

function [cliqcounts, sepcounts] = computeCounts(cliques, separators, npops, PARTITION)

ncliq = size(cliques,1);
nsep = size(separators,1);

cliqPartition = zeros(ncliq, size(cliques,2));
sepPartition = zeros(nsep, size(separators, 2));

apuCliq = find(cliques > 0);
apuSep = find(separators > 0);

cliqPartition(apuCliq) = PARTITION(cliques(apuCliq));
sepPartition(apuSep) = PARTITION(separators(apuSep));


cliqcounts = zeros(ncliq, npops);
for i = 1:npops
    cliqcounts(:,i) = sum(cliqPartition == i, 2);
end
    

sepcounts = zeros(nsep, npops);
for i = 1:npops
    sepcounts(:,i) = sum(sepPartition == i, 2);
end

%-------------------------------------------------------------------------

function diffInCliqCounts = computeDiffInCliqCounts(cliques, inds)
% Laskee muutoksen CLIQCOUNTS:ssa (tai SEPCOUNTS:ssa, jos syˆtteen?
% separators) kun yksilˆt inds siirret‰‰n.
% diffInCliqcounts on ncliq*1 taulu, joka on CLIQCOUNTS:n sarakkeesta josta
% yksilˆt inds siirret‰‰n ja lis‰tt‰v?sarakkeeseen, johon yksilˆt
% siirret‰‰n.

ncliq = size(cliques,1);
diffInCliqCounts = zeros(ncliq,1);
ninds = length(inds);
for i = 1:ninds
    ind = inds(i);
    rivit = sum((cliques == ind),2);
    diffInCliqCounts = diffInCliqCounts + rivit;
end

%-----------------------------------------------------------------------

function [logml, spatialPrior] = computeLogml(adjprior,priorTerm, ...
                                              CLIQCOUNTS, SEPCOUNTS, ...
                                              COUNTS, SUMCOUNTS)

notEmpty = any(CLIQCOUNTS);
npops = length(find(notEmpty == 1));
sumcliq=sum(CLIQCOUNTS, 2);
sumsep=sum(SEPCOUNTS, 2);
ncliq = size(CLIQCOUNTS, 1);
nsep = size(SEPCOUNTS, 1);

cliqsizes = sum(CLIQCOUNTS, 2)';
sepsizes = sum(SEPCOUNTS, 2)';
cliqsizes = min([cliqsizes; npops*ones(1,ncliq)])';
sepsizes = min([sepsizes; npops*ones(1,nsep)])';

klikkitn = sum(sum(gammaln(CLIQCOUNTS(:,notEmpty) + repmat(1./cliqsizes, [1 npops])))) ... 
                    - sum(npops*(gammaln(1./cliqsizes))) ...
                    - sum(gammaln(sumcliq + 1));
                                    
septn = sum(sum(gammaln(SEPCOUNTS(:,notEmpty) + repmat(1./sepsizes, [1 npops])))) ...
                - sum(npops*(gammaln(1./sepsizes))) ...
                - sum(gammaln(sumsep + 1));
            
                
%klikkitn = sum(sum(gammaln(CLIQCOUNTS + 1/npops))) ...
%                    - ncliq*npops*(gammaln(1/npops)) ...
%                    - sum(gammaln(sumcliq + 1));
%septn = sum(sum(gammaln(SEPCOUNTS + 1/npops))) ...
%                - nsep*npops*(gammaln(1/npops)) ...
%                - sum(gammaln(sumsep + 1));
            
spatialPrior = (klikkitn - septn);

%if spatialPrior > 0
%    keyboard
%end

x = size(COUNTS,1);
y = size(COUNTS,2);
z = size(COUNTS,3);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 z]) + COUNTS) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS),2) - priorTerm;    
    
logml = sum(popLogml) + spatialPrior;    
%--------------------------------------------------------------------------


function initializeGammaln(ninds, rowsFromInd, maxSize)
%Alustaa GAMMALN muuttujan s.e. GAMMALN(i,j)=gammaln((i-1) + 1/j)
global GAMMA_LN;
GAMMA_LN = zeros((1+ninds)*rowsFromInd, maxSize);
for i=1:(ninds+1)*rowsFromInd
    for j=1:maxSize
        GAMMA_LN(i,j)=gammaln((i-1) + 1/j);
    end
end


%----------------------------------------------------------------------------


function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sis‰lt‰‰ yksilˆiden inds2
% v‰liset et‰isyydet. ninds=kaikkien yksilˆiden lukum‰‰r?

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
rivi = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(rivi, 1) = inds2(i);
        apu(rivi, 2) = inds2(j);
        rivi = rivi+1;
    end
end
apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = dist(apu);


%----------------------------------------------------------------------------



function kunnossa = testaaGenePopData(tiedostonNimi)
% kunnossa == 0, jos data ei ole kelvollinen genePop data.
% Muussa tapauksessa kunnossa == 1.

kunnossa = 0;
fid = fopen(tiedostonNimi);
line1 = fgetl(fid);  %ensimm‰inen rivi
line2 = fgetl(fid);  %toinen rivi
line3 = fgetl(fid);  %kolmas

if (isequal(line1,-1) | isequal(line2,-1) | isequal(line3,-1))
    disp('Incorrect file format 1168'); fclose(fid);
    return 
end
if (testaaPop(line1)==1 || testaaPop(line2)==1)
    disp('Incorrect file format 1172'); fclose(fid);
    return
end
if testaaPop(line3)==1
    %2 rivi t‰llˆin lokusrivi
    nloci = rivinSisaltamienMjonojenLkm(line2);
    line4 = fgetl(fid);
    if isequal(line4,-1)
        disp('Incorrect file format 1180'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin nelj?t‰ytyy sis‰lt‰‰ pilkku.
        disp('Incorrect file format 1185'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedet‰‰n, ett?pys‰htyy
        pointer = pointer+1;
    end
    line4 = line4(pointer+1:end);  %pilkun j‰lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format 1195'); fclose(fid);
        return
    end
else
    line = fgetl(fid);
    lineNumb = 4;
    while (testaaPop(line)~=1 & ~isequal(line,-1))
        line = fgetl(fid);
        lineNumb = lineNumb+1;    
    end
    if isequal(line,-1)
        disp('Incorrect file format 1206'); fclose(fid);
        return
    end
    nloci = lineNumb-2;
    line4 = fgetl(fid);  %Eka rivi pop sanan j‰lkeen
    if isequal(line4,-1)
        disp('Incorrect file format 1212'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin t‰ytyy sis‰lt‰‰ pilkku.
        disp('Incorrect file format 1217'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedet‰‰n, ett?pys‰htyy.
        pointer = pointer+1;
    end
 
    line4 = line4(pointer+1:end);  %pilkun j‰lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format 1228'); fclose(fid);
        return
    end
end
kunnossa = 1;
fclose(fid);

%------------------------------------------------------
    

function [data, popnames] = lueGenePopData(tiedostonNimi)

fid = fopen(tiedostonNimi);
line = fgetl(fid);  %ensimm‰inen rivi
line = fgetl(fid);  %toinen rivi
count = rivinSisaltamienMjonojenLkm(line);

line = fgetl(fid);
lokusRiveja = 1;
while (testaaPop(line)==0)
    lokusRiveja = lokusRiveja+1;
    line = fgetl(fid);
end

if lokusRiveja>1
    nloci = lokusRiveja;
else
    nloci = count;
end

popnames = cell(10,2);
data = zeros(100, nloci+1);
nimienLkm=0;
ninds=0;
poimiNimi=1;
digitFormat = -1;
while line ~= -1
    line = fgetl(fid);
    
    if poimiNimi==1  
        %Edellinen rivi oli 'pop'
        nimienLkm = nimienLkm+1;
        ninds = ninds+1;
        if nimienLkm>size(popnames,1);
            popnames = [popnames; cell(10,2)];
        end
        nimi = lueNimi(line);
        if digitFormat == -1
            digitFormat = selvitaDigitFormat(line);
            divider = 10^digitFormat;
        end
        popnames{nimienLkm, 1} = {nimi};   %N‰in se on greedyMix:iss‰kin?!?
        popnames{nimienLkm, 2} = ninds;
        poimiNimi=0;
        
        data = addAlleles(data, ninds, line, divider);
        
    elseif testaaPop(line)
        poimiNimi = 1;
        
    elseif line ~= -1
        ninds = ninds+1;
        data = addAlleles(data, ninds, line, divider);
    end    
end

data = data(1:ninds*2,:);
popnames = popnames(1:nimienLkm,:);
fclose(fid);

npops = size(popnames,1);
ind = 1;
for pop = 1:npops
    if pop<npops
        while ind<popnames{pop+1,2}
            data([ind*2-1 ind*2],end) = pop;
            ind = ind+1;    
        end
    else
        while ind<=ninds
            data([ind*2-1 ind*2],end) = pop;
            ind = ind+1;    
        end
    end
end


%-------------------------------------------------------

function nimi = lueNimi(line)
%Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
n = 1;
merkki = line(n);
nimi = '';
while ~isequal(merkki,',')
    nimi = [nimi merkki];
    n = n+1;
    merkki = line(n);
end

%-------------------------------------------------------

function df = selvitaDigitFormat(line)
% line on ensimm‰inen pop-sanan j‰lkeinen rivi
% Genepop-formaatissa olevasta datasta. funktio selvitt‰‰
% rivin muodon perusteella, ovatko datan alleelit annettu
% 2 vai 3 numeron avulla.

n = 1;
merkki = line(n);
while ~isequal(merkki,',')
    n = n+1;
    merkki = line(n);
end

while ~any(merkki == '0123456789');
    n = n+1;
    merkki = line(n);
end
numeroja = 0;
while any(merkki == '0123456789');
    numeroja = numeroja+1;
    n = n+1;
    merkki = line(n);
end

df = numeroja/2;


%------------------------------------------------------


function count = rivinSisaltamienMjonojenLkm(line)
% Palauttaa line:n sis‰lt‰mien mjonojen lukum‰‰r‰n.
% Mjonojen v‰liss?t‰ytyy olla v‰lilyˆnti.
count = 0;
pit = length(line);
tila = 0;    %0, jos odotetaan v‰lilyˆntej? 1 jos odotetaan muita merkkej?
for i=1:pit
    merkki = line(i);
    if (isspace(merkki) & tila==0) 
        %Ei tehd?mit‰‰n.
    elseif (isspace(merkki) & tila==1)
        tila = 0;
    elseif (~isspace(merkki) & tila==0)
        tila = 1;
        count = count+1;
    elseif (~isspace(merkki) & tila==1)
        %Ei tehd?mit‰‰n
    end
end

%-------------------------------------------------------

function pal = testaaPop(rivi)
% pal=1, mik‰li rivi alkaa jollain seuraavista
% kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
% tapauksissa pal=0.

if length(rivi)<3 
    pal = 0;
    return
end
if (all(rivi(1:3)=='Pop') | ...
    all(rivi(1:3)=='pop') | ...
    all(rivi(1:3)=='POP'))
    pal = 1;    
    return
else 
    pal = 0;
    return
end


%--------------------------------------------------------


function data = addAlleles(data, ind, line, divider)
% Lisaa BAPS-formaatissa olevaan datataulukkoon
% yksilˆ‰ ind vastaavat rivit. Yksilˆn alleelit
% luetaan genepop-formaatissa olevasta rivist?
% line. Jos data on 3 digit formaatissa on divider=1000.
% Jos data on 2 digit formaatissa on divider=100.

nloci = size(data,2)-1;
if size(data,1) < 2*ind
    data = [data; zeros(100,nloci+1)];
end

k=1;
merkki=line(k);
while ~isequal(merkki,',')
   k=k+1;
   merkki=line(k);
end
line = line(k+1:end);
clear k; clear merkki;

alleeliTaulu = sscanf(line,'%d');

if length(alleeliTaulu)~=nloci
    disp('Incorrect data format.');
end

for j=1:nloci
    ekaAlleeli = floor(alleeliTaulu(j)/divider);
    if ekaAlleeli==0 ekaAlleeli=-999; end;
    tokaAlleeli = rem(alleeliTaulu(j),divider);
    if tokaAlleeli==0 tokaAlleeli=-999; end
    
    data(2*ind-1,j) = ekaAlleeli;
    data(2*ind,j) = tokaAlleeli;
end

data(2*ind-1,end) = ind;
data(2*ind,end) = ind;

%-------------------------------------------------------------------


function [varmuus,changesInLogml] = writeMixtureInfo(logml, globalRows, data, adjprior, ...
    priorTerm, outPutFile, inputFile, partitionSummary, popnames, ...
    cliques, separators, fixedK)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global LOGDIFF;

OUTPUT_FILE = 'baps5_output.baps';

ninds = size(globalRows,1);
npops =  size(COUNTS,3);
names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksilˆihin

if length(outPutFile)>0
    fid = fopen(outPutFile,'a');
else
    fid = -1;
    diary(OUTPUT_FILE); % save in text anyway.
end

dispLine;
disp('RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:');
disp(['Data file: ' inputFile]);
disp(['Number of clustered groups: ' ownNum2Str(ninds)]);
disp(['Number of clusters in optimal partition: ' ownNum2Str(npops)]);
disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
    fprintf(fid,'%s \n', ['RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clustered groups: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clusters in optimal partition: ' ownNum2Str(npops)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
end

cluster_count = length(unique(PARTITION));
disp(['Best Partition: ']);
if (fid ~= -1)
    fprintf(fid,'%s \n',['Best Partition: ']); fprintf(fid,'\n');
end
for m=1:cluster_count
    indsInM = find(PARTITION==m);
    length_of_beginning = 11 + floor(log10(m));
    cluster_size = length(indsInM);
    
    if names
        text = ['Cluster ' num2str(m) ': {' char(popnames{indsInM(1)})];
        for k = 2:cluster_size
            text = [text ', ' char(popnames{indsInM(k)})];
        end;    
    else
        text = ['Cluster ' num2str(m) ': {' num2str(indsInM(1))];
        for k = 2:cluster_size
            text = [text ', ' num2str(indsInM(k))];
        end;
    end
    text = [text '}'];
    while length(text)>58
        %Take one line and display it.
        new_line = takeLine(text,58);
        text = text(length(new_line)+1:end);
        disp(new_line);
        if (fid ~= -1)
            fprintf(fid,'%s \n',[new_line]);    
            fprintf(fid,'\n');
        end
        if length(text)>0
            text = [blanks(length_of_beginning) text];
        else
            text = [];
        end;
    end;
    if ~isempty(text)
        disp(text);
        if (fid ~= -1)
            fprintf(fid,'%s \n',[text]);
            fprintf(fid,'\n');
        end
    end
end

disp(' ');
disp(' ');
disp('Changes in log(marginal likelihood) if group i is moved to cluster j:');
if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if group i is moved to cluster j:']); fprintf(fid, '\n');
end

if names
    nameSizes = zeros(ninds,1);
    for i = 1:ninds
        nimi = char(popnames{i});
        nameSizes(i) = length(nimi);
    end
    maxSize = max(nameSizes);
    maxSize = max(maxSize, 5);
    erotus = maxSize - 5;
    alku = blanks(erotus);
    ekarivi = [alku 'group' blanks(6+erotus)];
else
    ekarivi = 'group      ';
end

for i = 1:cluster_count
    ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
end
disp(ekarivi);
if (fid ~= -1)
    fprintf(fid, '%s \n', [ekarivi]); fprintf(fid, '\n');
end

varmuus = zeros(ninds,1);
changesInLogml = LOGDIFF';

for ind = 1:ninds
    %[muutokset, diffInCounts] = laskeMuutokset(ind, globalRows, data, ...
    %    adjprior, priorTerm, logml, cliques, separators);
    %changesInLogml(:,ind) = muutokset;
    muutokset = changesInLogml(:,ind);
    if sum(exp(muutokset))>0
        varmuus(ind) = 1 - 1/sum(exp(muutokset));
    else
        varmuus(ind) = 0;
    end
    if names
        nimi = char(popnames{ind});
        rivi = [blanks(maxSize - length(nimi)) nimi ':'];
    else
        rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
    end
    for j = 1:npops
        rivi = [rivi '  ' logml2String(omaRound(muutokset(j)))];
    end
    disp(rivi);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
    end
end

disp(' '); disp(' ');
disp('KL-divergence matrix:');
dist_mat = zeros(npops, npops);
if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['KL-divergence matrix in PHYLIP format:']); %fprintf(fid, '\n');
end

maxnoalle = size(COUNTS,1);
nloci = size(COUNTS,2);
d = zeros(maxnoalle, nloci, npops);
prior = adjprior;
prior(find(prior==1))=0;
nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
prior(1,nollia)=1;
for pop1 = 1:npops
    d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
    %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
end
% ekarivi = blanks(7);
% for pop = 1:npops
%     ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
% end
ekarivi = num2str(npops);
disp(ekarivi);
if (fid ~= -1)
    fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
end

for pop1 = 1:npops
    rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
    for pop2 = 1:pop1-1
        dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
        div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
        div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
        div = (div12+div21)/2;
        % rivi = [rivi kldiv2str(div) '  '];
        dist_mat(pop1,pop2) = div;
    end
%     disp(rivi);
%     if (fid ~= -1)
%         fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
%     end
end

dist_mat = dist_mat + dist_mat'; % make it symmetric
for pop1 = 1:npops
    rivi = ['Cluster_' num2str(pop1) ' '];
    for pop2 = 1:npops
        rivi = [rivi kldiv2str(dist_mat(pop1,pop2)) ' '];
    end
    disp(rivi);
    if (fid ~= -1)
         fprintf(fid, '%s \n', [rivi]); %fprintf(fid, '\n');
    end
end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values']); fprintf(fid, '\n');
end

partitionSummaryKaikki = partitionSummary;
partitionSummary =[];
for i=1:size(partitionSummaryKaikki,3)
    partitionSummary = [partitionSummary; partitionSummaryKaikki(:,:,i)];
end
[I,J] = find(partitionSummaryKaikki(:,2,:)>-1e49);
partitionSummaryKaikki = partitionSummaryKaikki(I,:,:);
%keyboard
    
    
partitionSummary = sortrows(partitionSummary,2);
partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
partitionSummary = partitionSummary(find(partitionSummary(:,2)>-1e49),:);
if size(partitionSummary,1)>10
    vikaPartitio = 10;
else
    vikaPartitio = size(partitionSummary,1);
end
for part = 1:vikaPartitio
    line = [num2str(partitionSummary(part,1)) '    ' num2str(partitionSummary(part,2))];
    disp(line);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
    end
end

if ~fixedK

    disp(' ');
    disp(' ');
    disp('Probabilities for number of clusters');

    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['Probabilities for number of clusters']); fprintf(fid, '\n');
    end

    npopsTaulu = unique(partitionSummary(:,1));
    len = length(npopsTaulu);
    probs = zeros(len,1);
    partitionSummary(:,2) = partitionSummary(:,2)-max(partitionSummary(:,2));
    sumtn = sum(exp(partitionSummary(:,2)));
    for i=1:len
        npopstn = sum(exp(partitionSummary(find(partitionSummary(:,1)==npopsTaulu(i)),2)));
        probs(i) = npopstn / sumtn;
    end
    for i=1:len
        if probs(i)>1e-5
            line = [num2str(npopsTaulu(i)) '   ' num2str(probs(i))];
            disp(line);
            if (fid ~= -1)
                fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
            end
        end
    end
end

if (fid ~= -1)
    fclose(fid);
else
    diary off
end

%---------------------------------------------------------------


function dispLine;
disp('---------------------------------------------------');

%--------------------------------------------------------------------


function newline = takeLine(description,width)
%Returns one line from the description: line ends to the first
%space after width:th mark.
newLine = description(1:width);
n = width+1;
while ~isspace(description(n)) & n<length(description)
    n = n+1;
end;
newline = description(1:n);


%--------------------------------------------------------------

function num2 = omaRound(num)
% Pyˆrist‰‰ luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------


function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in?
% yks t‰ytyy olla kokonaisluku, joka on 
% v‰hint‰‰n -1:n suuruinen. Pienemmill?
% luvuilla tapahtuu jokin pyˆristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);


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


%--------------------------------------------------------------------------


function ninds = testaaOnkoKunnollinenBapsData(data)
%Tarkastaa onko viimeisess?sarakkeessa kaikki
%luvut 1,2,...,n johonkin n:‰‰n asti.
%Tarkastaa lis‰ksi, ett?on v‰hint‰‰n 2 saraketta.
if size(data,1)<2
    ninds = 0; return;
end
lastCol = data(:,end);
ninds = max(lastCol);
if ~isequal((1:ninds)',unique(lastCol))
    ninds = 0; return;
end


%--------------------------------------------------------------------------


function [ninds, data, heds] = testFastaData(inFile)
% added by Lu Cheng, 05.12.2012
if ~exist(inFile,'file')
    error('Fasta file %s does not exist!\n',inFile);
end

[heds, seqs]=fastaread(inFile);
ninds = length(seqs);

data = cell2mat(seqs(:));
newData = ones(size(data))*-9;
newData(ismember(data,'Aa'))=1;
newData(ismember(data,'Cc'))=2;
newData(ismember(data,'Gg'))=3;
newData(ismember(data,'Tt'))=4;
data = [newData (1:ninds)'];


%--------------------------------------------------------------------------

function [cliques, separators, vorPoints, vorCells, pointers] ...
    = handleCoords(coordinates)
%Laskee yksilˆiden luonnolliset naapurit koordinaateista.
%Naapurit lasketaan lis‰‰m‰ll?koordinaatteihin pisteit?
%jotta kutakin yksilˆ‰ vastaisi rajoitettu voronoi-solu
%Puuttuvat koordinaatit (negatiiviset) tulevat erakkopisteiksi
%
%M‰‰ritt‰‰ lis‰ksi yksilˆit?vastaavat voronoi tesselaation solut.
%vorPoints:ssa on solujen kulmapisteet ja vorCells:ss?kunkin solun
%kulmapisteiden indeksit. Pointers{i} sis‰lt‰‰ solussa i olevien yksilˆiden
%indeksit.



ninds = length(coordinates);
[I,J] = find(coordinates>0 | coordinates <0);  %K‰sitell‰‰n vain yksilˆit? joilta koordinaatit
I = unique(I);                %olemassa
ncoords = length(I);
puuttuvat = setdiff(1:ninds, I);
new_coordinates = addPoints(coordinates(I,:)); %Ymp‰rˆid‰‰n yksilˆt apupisteill?


apuData = [new_coordinates(1:ncoords,:) (1:ncoords)'];
apuData = sortrows(apuData,[1 2]);
erot = [diff(apuData(:,1)) diff(apuData(:,2))];
empties = find(erot(:,1)==0 & erot(:,2)==0);
samat = cell(length(empties),1);
pointer = 0;

for i = 1:length(empties)
    if i == 1 | empties(i) - empties(i-1) > 1  %Tutkitaan onko eri pisteess?kuin edellinen
        pointer = pointer+1;
        samat{pointer} = [apuData(empties(i),3) apuData(empties(i)+1,3)];
    else
        samat{pointer} = [samat{pointer} apuData(empties(i)+1,3)];
    end
end

samat = samat(1:pointer);

erot = []; apuData = []; empties = [];

%tri = delaunay(new_coordinates(:,1), new_coordinates(:,2), {'Qt','Qbb','Qc','Qz'});    %Apupisteiden takia ok.
tri = delaunay(new_coordinates(:,1), new_coordinates(:,2)); 
%[rivi,sarake] = find(tri>ncoords);    %J‰tet‰‰n huomiotta apupisteet
%tri(rivi,:) = [];
pituus = tri(:,1);
pituus = length(pituus);
parit = zeros(6*pituus,2);
for i = 1:pituus                        %Muodostetaan kolmikoista parit
    j = 6*(i-1)+1;
    parit(j,:) = tri(i,1:2);
    parit(j+1,:) = tri(i,1:2:3);
    parit(j+2,:) = tri(i,2:3);
    parit(j+3:j+5,:) = [parit(j:j+2,2) parit(j:j+2,1)];
end
parit = unique(parit,'rows');
[rivi,sarake] = find(parit>ncoords);     %J‰tet‰‰n huomiotta apupisteet
parit(rivi,:) = [];
parit = I(parit);                         %Otetaan poistetut takaisin mukaan
graph = sparse(parit(:,1),parit(:,2),1, ninds, ninds);


%Kopioidaan samassa pisteess?olevien yksilˆiden naapurustot
%silt? jolle ne laitettu.

    for i = 1:length(samat);
        taulu = I(samat{i});        
        [rivi,sarake] = find(graph(taulu,:)>0);
        if length(rivi) > 0
            kopioitava = graph(taulu(rivi(1)),:);
            for j = 1:length(taulu);
                graph(taulu(j),:) = kopioitava;
                graph(:,taulu(j)) = kopioitava';
            end
        end
    end

    %Asetetaan samassa pisteess?olevat yksilˆt toistensa naapureiksi

    for i = 1:length(samat)
        for j = I(samat{i})
            for k = I(samat{i})
                if k ~= j
                    graph(j,k) = 1;
                end
            end
        end
    end

%Laskee maksimin klikkien ja separaattorien koolle
%M‰‰ritet‰‰n myˆs klikit ja separaattorit 

[ncliq, nsep, cliq, sep] = laskeKlikit(graph, ninds, ninds);

sumcliq = sum(ncliq);
sumsep = sum(nsep);
maxCliqSize = max(find(sumcliq > 0));
maxSepSize = max(find(sumsep > 0));

cliques = zeros(length(cliq), maxCliqSize);
separators = zeros(length(sep), maxSepSize);

nollia = zeros(1, length(cliq));
for i = 1:length(cliq);
    klikki = cliq{i};
    if length(klikki)>1
        cliques(i, 1:length(klikki)) = klikki;
    else
        nollia(i)=1;
    end
end
cliques(find(nollia==1), :) = [];

for i = 1:length(sep);
    klikki = sep{i};
    separators(i, 1:length(klikki)) = klikki;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M‰‰ritet‰‰n yksilˆit?vastaavat voronoi tesselaation solut

[vorPoints, vorCells] = voronoin(new_coordinates, {'Qbb', 'Qz'});

bounded = ones(length(vorCells),1);
for i=1:length(vorCells)
    if isempty(vorCells{i}) || length(find(vorCells{i}==1))>0
        bounded(i)=0;
    end
end



vorCells = vorCells(find(bounded == 1));
pointers = cell(length(vorCells),1);
empties = zeros(1,length(vorCells));
X = coordinates(:,1);
Y = coordinates(:,2);

for i=1:length(pointers)
    vx = vorPoints(vorCells{i},1);
    vy = vorPoints(vorCells{i},2);
    IN = inpolygon(X,Y,vx,vy);
    if any(IN)==0
        empties(i) = 1;
    else
        pointers{i} = find(IN ==1)';
    end
end

vorCells = vorCells(find(empties == 0));
pointers = pointers(find(empties == 0));

%--------------------------------------------------------------------------

function [ncliques, nseparators, cliques, separators] = ...
    laskeKlikit(M, maxCliqSize,maxSepSize)
%Laskee samankokoisten klikkien m‰‰r‰n verkosta M
%ncliques(i)=kokoa i olevien klikkien m‰‰r?
%nseparators vastaavasti

ncliques=zeros(1,maxCliqSize); 
nseparators=zeros(1,maxSepSize);

if isequal(M,[])
    return;
end

[cliques,separators]=findCliques(M);

for i=1:length(cliques)
    ncliques(length(cliques{i}))=ncliques(length(cliques{i}))+1;
end

%cliqmax=max(find(ncliques~=0));
%ncliques=ncliques(1:cliqmax);

for i=1:length(separators)
    nseparators(length(separators{i}))=nseparators(length(separators{i}))+1;
end

%sepmax=max(find(nseparators~=0));
%nseparators=nseparators(1:sepmax);

%--------------------------------------------------------------------------

function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }
%
% Original by Kevin Murphy, modified by Leon Peshkin

if isempty(A)
    C = [];
    return;
elseif isempty(B)
    C = A;
    return; 
else % both non-empty
    bits = zeros(1, max(max(A), max(B)));
    bits(A) = 1;
    bits(B) = 0;
    C = A(logical(bits(A)));
end


%--------------------------------------------------------------------------

function logml = checkLogml(priorTerm, adjprior, cliques, separators)
% tarkistaa logml:n

global CLIQCOUNTS;
global SEPCOUNTS;
global PARTITION;

npops = length(unique(PARTITION));
[cliqcounts, sepcounts] = computeCounts(cliques, separators, npops);
        
CLIQCOUNTS = cliqcounts;
SEPCOUNTS = sepcounts;
    

[logml, spatialPrior] = computeLogml(adjprior, priorTerm);
    
disp(['logml: ' logml2String(logml) ', spatial prior: ' logml2String(spatialPrior)]);

%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimm‰isen tyhj‰n populaation indeksin. Jos tyhji?
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;

pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end

%--------------------------------------------------------------------------

% function viallinen = testaaKoordinaatit(ninds, coordinates)
% % Testaa onko koordinaatit kunnollisia.
% 
% viallinen = 1;
% if ~isnumeric(coordinates)
%     return
% end
% 
% oikeanKokoinen = (size(coordinates,1) == ninds) & (size(coordinates,2) == 2);
% if oikeanKokoinen
%     viallinen = 0;
% end

function [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates)
% Testaa onko koordinaatit kunnollisia.
% modified by Lu Cheng, 05.12.2012

viallinen = 1;
if ~isnumeric(coordinates)
    warning('Coordinates are not numerical!');
    return;
end

oikeanKokoinen = (size(coordinates,1) == ninds) & (size(coordinates,2) == 2);
if ~oikeanKokoinen
    warning('Wrong coordinates dimension!');
    return;
end

posstr = cellfun(@(x) sprintf('%.10f',x),num2cell(coordinates),'UniformOutput',false);
posstr = cellfun(@(x) regexprep(x,'0+$',''),posstr,'UniformOutput',false);

uni1 = unique(posstr(:,1));
uni2 = unique(posstr(:,2));
posstr_new = posstr;

if length(uni1)==ninds && length(uni2)==ninds
    viallinen = 0;
    return;
else
    ans = questdlg('Input coordinates are not unique. Do you want to make them unique?','coordinates NOT unique', 'Yes','No','Yes');
    if strcmp(ans,'No')
        warning('Coordinates are not unique!');
        return;
    end 
end

for i=1:length(uni1)
    tmpinds = find(ismember(posstr(:,1),uni1(i)));
    tmpNinds = length(tmpinds);
    
    if tmpNinds==1
        continue;
    end
    
    assert(tmpNinds<100);
    tmparr = round(linspace(0,99,tmpNinds));
    tmparr = tmparr(randperm(tmpNinds));
    
    for j=1:tmpNinds
        posstr_new{tmpinds(j),1}=sprintf('%s%02d',posstr{tmpinds(j),1},tmparr(j));
    end
end

for i=1:length(uni2)
    tmpinds = find(ismember(posstr(:,2),uni2(i)));
    tmpNinds = length(tmpinds);
    
    if tmpNinds==1
        continue;
    end
    
    assert(tmpNinds<100);
    tmparr = round(linspace(0,99,tmpNinds));
    tmparr = tmparr(randperm(tmpNinds));
    
    for j=1:tmpNinds
        posstr_new{tmpinds(j),2}=sprintf('%s%02d',posstr{tmpinds(j),2},tmparr(j));
    end
end

coordinates = cellfun(@str2double,posstr_new);
uni1 = unique(coordinates(:,1));
uni2 = unique(coordinates(:,2));
if length(uni1)==ninds && length(uni2)==ninds
    viallinen = 0;
else
    warning('Can not make coordinates unique!');
end

%--------------------------------------------------------------------------


function [sumcounts, counts] = ...
    initialCounts(partition, data, npops, noalle)

nloci=size(data,2);
% ninds = size(rows, 1);

%koot = rows(:,1) - rows(:,2) + 1;
%maxSize = max(koot);

counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);
for i=1:npops
    for j=1:nloci
        havainnotLokuksessa = find(partition==i & data(:,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end


%--------------------------------------------------------------------------

function [popnames2, rowsFromInd] = findOutRowsFromInd(popnames, rows)

ploidisuus = questdlg('Specify the type of individuals in the data: ',...
    'Individual type?', 'Haploid', 'Diploid', 'Tetraploid', ...
    'Diploid');

switch ploidisuus    
case 'Haploid'
    rowsFromInd = 1;
case 'Diploid'
    rowsFromInd = 2;
case 'Tetraploid'
    rowsFromInd = 4;
end

if ~isempty(popnames)
    for i = 1:size(rows,1)
        popnames2{i,1} = popnames{i,1};
        rivi = rows(i,1):rows(i,2);
        popnames2{i,2} = (rivi(rowsFromInd))/rowsFromInd;    
    end
else 
    popnames2 = [];
end

%--------------------------------------------------------------------------

function fiksaaPartitioYksiloTasolle(rows, rowsFromInd)

global PARTITION;
totalRows = 0;
for ind = 1:size(rows,1)
    totalRows = totalRows + (rows(ind,2)-rows(ind,1)+1);
end
partitio2 = zeros(totalRows/rowsFromInd,1);

for ind = 1:size(rows,1)
    kaikkiRivit = rows(ind,1):rows(ind,2);
    for riviNumero = rowsFromInd:rowsFromInd:length(kaikkiRivit)
    %for riviNumero = rowsFromInd:rowsFromInd:length(rows{ind})
        %rivi = rows{ind}(riviNumero);
        rivi = kaikkiRivit(riviNumero);
        partitio2(rivi/rowsFromInd) = PARTITION(ind);
    end
end
PARTITION = partitio2;
