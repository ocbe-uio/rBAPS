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
    a_data = data(:,1:end-1);

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

%-------------------------------------------------------------------


function changesInLogml = writeMixtureInfoPop(logml, rows, data, adjprior, ...
    priorTerm, outPutFile, inputFile, partitionSummary, popnames, fixedK)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global LOGDIFF;
ninds = size(rows,1);
npops =  size(COUNTS,3);
names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksilöihin
changesInLogml = [];
if length(outPutFile)>0
    fid = fopen(outPutFile,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
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
    fprintf(fid,'\n');
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
    end;
end

if npops > 1
    disp(' ');
    disp(' ');
    disp('Changes in log(marginal likelihood) if group i is moved to cluster j:');
    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if group i is moved to cluster j:']); %fprintf(fid, '\n');
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
        fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
    end

    changesInLogml = LOGDIFF';
    for ind = 1:ninds
        %[muutokset, diffInCounts] = laskeMuutokset(ind, rows, data, ...
        %    adjprior, priorTerm);
        %changesInLogml(:,ind) = muutokset;
        muutokset = changesInLogml(:,ind);
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
    disp('KL-divergence matrix in PHYLIP format:');
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
%     ekarivi = blanks(7);
%     for pop = 1:npops
%         ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
%     end
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
%         disp(rivi);
%         if (fid ~= -1)
%             fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
%         end
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

end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values']); fprintf(fid, '\n');
end

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
