function greedyPopMix

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
clearGlobalVars;

% check whether fixed k mode is selected
h0 = findobj('Tag','fixk_menu');
fixedK = get(h0, 'userdata');

if fixedK
    if ~(fixKWarning == 1) % call function fixKWarning
        return
    end
end

% check whether partition compare mode is selected
h1 = findobj('Tag','partitioncompare_menu');
partitionCompare = get(h1, 'userdata');

% LASKENNAN ALKUARVOJEN MÄÄRITTÄMINEN

input_type = questdlg('Specify the format of your data: ',...
    'Specify Data Format', ...
    'BAPS-format', 'GenePop-format', 'Preprocessed data', ...
    'BAPS-format');

if isempty(input_type)
    return
end

if isequal(input_type,'BAPS-format')  %Raakadata
    waitALittle;
    [filename, pathname] = uigetfile('*.txt', 'Load data in BAPS-format');
    if filename==0
        return;
    end
    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
    end
    data = load([pathname filename]);
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    if (ninds==0)
        disp('Incorrect Data-file.');
        return;
    end
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    rowsFromInd = 0;  %Ei tiedet?
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;

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

elseif isequal(input_type,'GenePop-format')
    waitALittle;
    [filename, pathname] = uigetfile('*.txt', 'Load data in GenePop-format');
    if filename==0
        return;
    end
    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
    end
    kunnossa = testaaGenePopData([pathname filename]);
    if kunnossa==0
        return
    end

    [data, popnames]=lueGenePopDataPop([pathname filename]);
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    rowsFromInd = 2;  %Tiedetään GenePop:in tapauksessa.

    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;
end

if ~isequal(input_type, 'Preprocessed data')
    2

    npops = size(rows,1);
    PARTITION = 1:npops';  %Jokainen "yksil? eli populaatio on oma ryhmäns?
    [sumcounts, counts, logml] = ...
        initialPopCounts(a_data, npops, rows, noalle, adjprior);
    COUNTS = counts; SUMCOUNTS = sumcounts;
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);

    clear('counts', 'sumcounts','pathname','filename','vast2',...
        'vast3','vast4');
    [Z,dist] = getPopDistancesByKL(adjprior);  %Saadaan COUNTS:in avulla.

    save_preproc = questdlg('Do you wish to save pre-processed data?',...
        'Save pre-processed data?',...
        'Yes','No','Yes');
    if isequal(save_preproc,'Yes');
        waitALittle;
        [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
        kokonimi = [pathname filename];
        c.data = data; c.rows = rows; c.alleleCodes = alleleCodes;
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.dist = dist; c.Z = Z; c.popnames = popnames; c.rowsFromInd = rowsFromInd;
        c.npops = npops;  c.logml = logml;
%         save(kokonimi,'c');
        save(kokonimi,'c','-v7.3'); % Lu Cheng, 08.06.2012
        clear c;
    end;
end

if isequal(input_type, 'Preprocessed data')
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
    %load([pathname filename],'c');
    %if ~exist('c')   %TESTAUS
    %    disp('Incorrect file format.');
    %    return
    %elseif ~isfield(c,'rows')
    %    disp('Incorrect file format.');
    %    return
    %end
    struct_array = load([pathname filename]);
    if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        if ~isfield(c,'rows')
            disp('Incorrect file format');
            return
        end
    elseif isfield(struct_array,'rows')  %Mideva versio
        c = struct_array;
    else
        disp('Incorrect file format');
        return;
    end
    data = double(c.data); rows = c.rows; alleleCodes = c.alleleCodes;
    noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
    dist = c.dist; Z = c.Z; popnames = c.popnames; rowsFromInd = c.rowsFromInd;
    clear c;
end

c.data=data; c.rows = rows; c.alleleCodes = alleleCodes;
c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
c.dist=dist; c.Z=Z; c.rowsFromInd = rowsFromInd;

% partition compare
if ~isempty(partitionCompare)
    nsamplingunits = size(rows,1);
    partitions = partitionCompare.partitions;
    npartitions = size(partitions,2);
    partitionLogml = zeros(1,npartitions);
    for i = 1:npartitions
        % number of unique partition lables
        npops = length(unique(partitions(:,i)));
        try
            partitionInd = zeros(rows(end),1);
            partitionSample = partitions(:,i);
            for j = 1: nsamplingunits
                partitionInd([c.rows(j,1):c.rows(j,2)]) = partitionSample(j);
            end
            partitionLogml(i) = ...
                initialCounts(partitionInd, data(:,1:end-1), npops, c.rows, noalle, adjprior);
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
    [logml, npops, partitionSummary]=indMix_fixK(c);
else
    [logml, npops, partitionSummary]=indMix(c);
end

if logml==1
    return
end

data = data(:,1:end-1);
viewPopMixPartition(PARTITION, rows, popnames);
%npops = poistaTyhjatPopulaatiot(npops);
%POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);

h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
h0 = findobj('Tag','filename2_text');
outp = get(h0,'String');
changesInLogml = writeMixtureInfoPop(logml, rows, data, adjprior, priorTerm, ...
    outp,inp,partitionSummary, popnames, fixedK);

talle = questdlg(['Do you want to save the mixture populations ' ...
    'so that you can use them later in admixture analysis?'], ...
    'Save results?','Yes','No','Yes');
if isequal(talle,'Yes')
    waitALittle;
    [filename, pathname] = uiputfile('*.mat','Save results as');

    if (filename == 0) & (pathname == 0)
        % Cancel was pressed
        return
    else % copy 'baps4_output.baps' into the text file with the same name.
        if exist('baps4_output.baps','file')
            copyfile('baps4_output.baps',[pathname filename '.txt'])
            delete('baps4_output.baps')
        end
    end

    if rowsFromInd==0
        %Käytettiin BAPS-formaattia, eik?rowsFromInd ole tunnettu.
        [popnames, rowsFromInd] = findOutRowsFromInd(popnames, rows);
    end

    groupPartition = PARTITION;

    fiksaaPartitioYksiloTasolle(rows, rowsFromInd);

    c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
    c.alleleCodes = alleleCodes; c.adjprior = adjprior;
    c.rowsFromInd = rowsFromInd; c.popnames = popnames;
    c.data = data; c.npops = npops; c.noalle = noalle;
    c.mixtureType = 'popMix'; c.groupPartition = groupPartition;
    c.rows = rows; c.logml = logml; c.changesInLogml = changesInLogml;
%     save([pathname filename], 'c');
    save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
else
    if exist('baps4_output.baps','file')
        delete('baps4_output.baps')
    end
end
