----------------------- line 1 -----------------------
function greedyPopMix
function greedyPopMix
----------------------- line 2 -----------------------


----------------------- line 3 -----------------------
global PARTITION;
# global PARTITION
----------------------- line 4 -----------------------
global COUNTS;
# global COUNTS
----------------------- line 5 -----------------------
global SUMCOUNTS;
# global SUMCOUNTS
----------------------- line 6 -----------------------
global POP_LOGML;
# global POP_LOGML
----------------------- line 7 -----------------------
clearGlobalVars;
clearGlobalVars
----------------------- line 8 -----------------------


----------------------- line 9 -----------------------
% check whether fixed k mode is selected
#  check whether fixed k mode is selected
----------------------- line 10 -----------------------
h0 = findobj('Tag','fixk_menu');
h0 <- findobj('Tag', 'fixk_menu')
----------------------- line 11 -----------------------
fixedK = get(h0, 'userdata');
fixedK <- get(h0, 'userdata')
----------------------- line 12 -----------------------


----------------------- line 13 -----------------------
if fixedK
if (fixedK) {
----------------------- line 14 -----------------------
    if ~(fixKWarning == 1) % call function fixKWarning
    if (!(fixKWarning == 1)# call function fixKWarning) {
----------------------- line 15 -----------------------
        return
        return
----------------------- line 16 -----------------------
    end
    }
----------------------- line 17 -----------------------
end
}
----------------------- line 18 -----------------------


----------------------- line 19 -----------------------
% check whether partition compare mode is selected
#  check whether partition compare mode is selected
----------------------- line 20 -----------------------
h1 = findobj('Tag','partitioncompare_menu');
h1 <- findobj('Tag', 'partitioncompare_menu')
----------------------- line 21 -----------------------
partitionCompare = get(h1, 'userdata');
partitionCompare <- get(h1, 'userdata')
----------------------- line 22 -----------------------


----------------------- line 23 -----------------------
% LASKENNAN ALKUARVOJEN M��RITT�MINEN
#  LASKENNAN ALKUARVOJEN M < c4 > < c4 > RITT < c4 > MINEN
----------------------- line 24 -----------------------


----------------------- line 25 -----------------------
input_type = questdlg('Specify the format of your data: ',...
input_type <- questdlg('Specify the format of your data: ',
----------------------- line 26 -----------------------
    'Specify Data Format', ...
    'Specify Data Format',
----------------------- line 27 -----------------------
    'BAPS-format', 'GenePop-format', 'Preprocessed data', ...
    'BAPS - format', 'GenePop - format', 'Preprocessed data',
----------------------- line 28 -----------------------
    'BAPS-format');
    'BAPS - format')
----------------------- line 29 -----------------------


----------------------- line 30 -----------------------
if isempty(input_type)
if (is.null(input_type)) {
----------------------- line 31 -----------------------
    return
    return
----------------------- line 32 -----------------------
end
}
----------------------- line 33 -----------------------


----------------------- line 34 -----------------------
if isequal(input_type,'BAPS-format')  %Raakadata
if (isequal(input_type, 'BAPS - format') # Raakadata) {
----------------------- line 35 -----------------------
    waitALittle;
    waitALittle
----------------------- line 36 -----------------------
    [filename, pathname] = uigetfile('*.txt', 'Load data in BAPS-format');
    [filename, pathname] = uigetfile(' * .txt', 'Load data in BAPS - format')
----------------------- line 37 -----------------------
    if filename==0
    if (filename == 0) {
----------------------- line 38 -----------------------
        return;
        return
----------------------- line 39 -----------------------
    end
    }
----------------------- line 40 -----------------------
    if ~isempty(partitionCompare)
    if (!is.null(partitionCompare)) {
----------------------- line 41 -----------------------
        fprintf(1,'Data: %s\n',[pathname filename]);
        fprintf(1, 'Data:# s\n', [pathname filename])
----------------------- line 42 -----------------------
    end
    }
----------------------- line 43 -----------------------
    data = load([pathname filename]);
    data <- load([pathname filename])
----------------------- line 44 -----------------------
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    ninds <- testaaOnkoKunnollinenBapsData(data) # TESTAUS
----------------------- line 45 -----------------------
    if (ninds==0)
    if ((ninds == 0)) {
----------------------- line 46 -----------------------
        disp('Incorrect Data-file.');
        disp('Incorrect Data - file.')
----------------------- line 47 -----------------------
        return;
        return
----------------------- line 48 -----------------------
    end
    }
----------------------- line 49 -----------------------
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data)
----------------------- line 50 -----------------------
    rowsFromInd = 0;  %Ei tiedet?
    rowsFromInd <- 0 # Ei tiedet?
----------------------- line 51 -----------------------
    h0 = findobj('Tag','filename1_text');
    h0 <- findobj('Tag', 'filename1_text')
----------------------- line 52 -----------------------
    set(h0,'String',filename); clear h0;
    set(h0, 'String', filename) clear h0
----------------------- line 53 -----------------------


----------------------- line 54 -----------------------
    load_names = questdlg('Do you wish to specify the names of the groups?',...
    load_names <- questdlg('Do you wish to specify the names of the groups?',
----------------------- line 55 -----------------------
        'Input group names?','Yes','No','Yes');
        'Input group names?', 'Yes', 'No', 'Yes')
----------------------- line 56 -----------------------
    if isequal(load_names,'Yes')
    if (isequal(load_names, 'Yes')) {
----------------------- line 57 -----------------------
        waitALittle;
        waitALittle
----------------------- line 58 -----------------------
        [filename, pathname] = uigetfile('*.txt', 'Load group names');
        [filename, pathname] = uigetfile(' * .txt', 'Load group names')
----------------------- line 59 -----------------------
        popnames = initPopNames([pathname filename]);
        popnames <- initPopNames([pathname filename])
----------------------- line 60 -----------------------
        if (size(popnames,1)~=ninds)
        if ((size[popnames, 1] <- ninds)) {
----------------------- line 61 -----------------------
            disp('Incorrect name-file.');
            disp('Incorrect name-file.')
----------------------- line 62 -----------------------
            popnames = [];
            popnames <- vector()
----------------------- line 63 -----------------------
        end
        }
----------------------- line 64 -----------------------
    else
    } else {
----------------------- line 65 -----------------------
        popnames = [];
        popnames <- vector()
----------------------- line 66 -----------------------
    end
    }
----------------------- line 67 -----------------------


----------------------- line 68 -----------------------
elseif isequal(input_type,'GenePop-format')
} else if (isequal(input_type, 'GenePop - format')) {
----------------------- line 69 -----------------------
    waitALittle;
    waitALittle
----------------------- line 70 -----------------------
    [filename, pathname] = uigetfile('*.txt', 'Load data in GenePop-format');
    [filename, pathname] = uigetfile(' * .txt', 'Load data in GenePop - format')
----------------------- line 71 -----------------------
    if filename==0
    if (filename == 0) {
----------------------- line 72 -----------------------
        return;
        return
----------------------- line 73 -----------------------
    end
    }
----------------------- line 74 -----------------------
    if ~isempty(partitionCompare)
    if (!is.null(partitionCompare)) {
----------------------- line 75 -----------------------
        fprintf(1,'Data: %s\n',[pathname filename]);
        fprintf(1, 'Data:# s\n', [pathname filename])
----------------------- line 76 -----------------------
    end
    }
----------------------- line 77 -----------------------
    kunnossa = testaaGenePopData([pathname filename]);
    kunnossa <- testaaGenePopData([pathname filename])
----------------------- line 78 -----------------------
    if kunnossa==0
    if (kunnossa == 0) {
----------------------- line 79 -----------------------
        return
        return
----------------------- line 80 -----------------------
    end
    }
----------------------- line 81 -----------------------


----------------------- line 82 -----------------------
    [data, popnames]=lueGenePopDataPop([pathname filename]);
    [data, popnames]=lueGenePopDataPop([pathname filename])
----------------------- line 83 -----------------------
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data)
----------------------- line 84 -----------------------
    rowsFromInd = 2;  %Tiedet��n GenePop:in tapauksessa.
    rowsFromInd <- 2 # Tiedet < e4 > < e4 > n GenePop:in tapauksessa.
----------------------- line 85 -----------------------


----------------------- line 86 -----------------------
    h0 = findobj('Tag','filename1_text');
    h0 <- findobj('Tag', 'filename1_text')
----------------------- line 87 -----------------------
    set(h0,'String',filename); clear h0;
    set(h0, 'String', filename) clear h0
----------------------- line 88 -----------------------
end
}
----------------------- line 89 -----------------------


----------------------- line 90 -----------------------
if ~isequal(input_type, 'Preprocessed data')
if (!isequal(input_type, 'Preprocessed data')) {
----------------------- line 91 -----------------------
    a_data = data(:,1:end-1);
    a_data <- data[, 1:end - 1]
----------------------- line 92 -----------------------


----------------------- line 93 -----------------------
    npops = size(rows,1);
    npops <- size(rows, 1)
----------------------- line 94 -----------------------
    PARTITION = 1:npops';  %Jokainen "yksil? eli populaatio on oma ryhm�ns?
    PARTITION <- 1:npops' # Jokainen "yksil? eli populaatio on oma ryhm < e4 > ns?
----------------------- line 95 -----------------------
    [sumcounts, counts, logml] = ...
    [sumcounts, counts, logml] =
----------------------- line 96 -----------------------
        initialPopCounts(a_data, npops, rows, noalle, adjprior);
        initialPopCounts(a_data, npops, rows, noalle, adjprior)
----------------------- line 97 -----------------------
    COUNTS = counts; SUMCOUNTS = sumcounts;
    COUNTS <- counts SUMCOUNTS <- sumcounts
----------------------- line 98 -----------------------
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    POP_LOGML <- computePopulationLogml(1:npops, adjprior, priorTerm)
----------------------- line 99 -----------------------


----------------------- line 100 -----------------------
    clear('counts', 'sumcounts','pathname','filename','vast2',...
    clear('counts', 'sumcounts', 'pathname', 'filename', 'vast2',
----------------------- line 101 -----------------------
        'vast3','vast4');
        'vast3', 'vast4')
----------------------- line 102 -----------------------
    [Z,dist] = getPopDistancesByKL(adjprior);  %Saadaan COUNTS:in avulla.
    [Z, dist] = getPopDistancesByKL(adjprior) # Saadaan COUNTS:in avulla.
----------------------- line 103 -----------------------


----------------------- line 104 -----------------------
    save_preproc = questdlg('Do you wish to save pre-processed data?',...
    save_preproc <- questdlg('Do you wish to save pre-processed data?',
----------------------- line 105 -----------------------
        'Save pre-processed data?',...
        'Save pre-processed data?',
----------------------- line 106 -----------------------
        'Yes','No','Yes');
        'Yes', 'No', 'Yes')
----------------------- line 107 -----------------------
    if isequal(save_preproc,'Yes');
    if (isequal(save_preproc, 'Yes')) {
----------------------- line 108 -----------------------
        waitALittle;
        waitALittle
----------------------- line 109 -----------------------
        [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
        [filename, pathname] = uiputfile(' * .mat', 'Save pre-processed data as')
----------------------- line 110 -----------------------
        kokonimi = [pathname filename];
        kokonimi = [pathname filename]
----------------------- line 111 -----------------------
        c.data = data; c.rows = rows; c.alleleCodes = alleleCodes;
        c.data <- data c.rows <- rows c.alleleCodes <- alleleCodes
----------------------- line 112 -----------------------
        c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
        c.noalle <- noalle c.adjprior <- adjprior c.priorTerm <- priorTerm
----------------------- line 113 -----------------------
        c.dist = dist; c.Z = Z; c.popnames = popnames; c.rowsFromInd = rowsFromInd;
        c.dist <- dist c.Z <- Z c.popnames <- popnames c.rowsFromInd <- rowsFromInd
----------------------- line 114 -----------------------
        c.npops = npops;  c.logml = logml;
        c.npops <- npops  c.logml <- logml
----------------------- line 115 -----------------------
%         save(kokonimi,'c');
#          save(kokonimi, 'c')
----------------------- line 116 -----------------------
        save(kokonimi,'c','-v7.3'); % Lu Cheng, 08.06.2012
        save(kokonimi, 'c', ' - v7.3')# Lu Cheng, 08.06.2012
----------------------- line 117 -----------------------
        clear c;
        clear c
----------------------- line 118 -----------------------
    end;
    }
----------------------- line 119 -----------------------
end
}
----------------------- line 120 -----------------------


----------------------- line 121 -----------------------
if isequal(input_type, 'Preprocessed data')
if (isequal(input_type, 'Preprocessed data')) {
----------------------- line 122 -----------------------
    waitALittle;
    waitALittle
----------------------- line 123 -----------------------
    [filename, pathname] = uigetfile('*.mat', 'Load pre-processed data');
    [filename, pathname] = uigetfile(' * .mat', 'Load pre-processed data')
----------------------- line 124 -----------------------
    if filename==0
    if (filename == 0) {
----------------------- line 125 -----------------------
        return;
        return
----------------------- line 126 -----------------------
    end
    }
----------------------- line 127 -----------------------


----------------------- line 128 -----------------------
    if ~isempty(partitionCompare)
    if (!is.null(partitionCompare)) {
----------------------- line 129 -----------------------
        fprintf(1,'Data: %s\n',[pathname filename]);
        fprintf(1, 'Data:# s\n', [pathname filename])
----------------------- line 130 -----------------------
    end
    }
----------------------- line 131 -----------------------


----------------------- line 132 -----------------------
    h0 = findobj('Tag','filename1_text');
    h0 <- findobj('Tag', 'filename1_text')
----------------------- line 133 -----------------------
    set(h0,'String',filename); clear h0;
    set(h0, 'String', filename) clear h0
----------------------- line 134 -----------------------
    %load([pathname filename],'c');
   # load([pathname filename], 'c')
----------------------- line 135 -----------------------
    %if ~exist('c')   %TESTAUS
   # if (!exist('c')  # TESTAUS) {
----------------------- line 136 -----------------------
    %    disp('Incorrect file format.');
#     disp('Incorrect file format.')
----------------------- line 137 -----------------------
    %    return
#     return
----------------------- line 138 -----------------------
    %elseif ~isfield(c,'rows')
# } else if (!isfield(c, 'rows')) {
----------------------- line 139 -----------------------
    %    disp('Incorrect file format.');
#     disp('Incorrect file format.')
----------------------- line 140 -----------------------
    %    return
#     return
----------------------- line 141 -----------------------
    %end
# }
----------------------- line 142 -----------------------
    struct_array = load([pathname filename]);
    struct_array <- load([pathname filename])
----------------------- line 143 -----------------------
    if isfield(struct_array,'c')  %Matlab versio
    if (isfield(struct_array, 'c') # Matlab versio) {
----------------------- line 144 -----------------------
        c = struct_array.c;
        c <- struct_array.c
----------------------- line 145 -----------------------
        if ~isfield(c,'rows')
        if (!isfield(c, 'rows')) {
----------------------- line 146 -----------------------
            disp('Incorrect file format');
            disp('Incorrect file format')
----------------------- line 147 -----------------------
            return
            return
----------------------- line 148 -----------------------
        end
        }
----------------------- line 149 -----------------------
    elseif isfield(struct_array,'rows')  %Mideva versio
    } else if (isfield(struct_array, 'rows') # Mideva versio) {
----------------------- line 150 -----------------------
        c = struct_array;
        c <- struct_array
----------------------- line 151 -----------------------
    else
    } else {
----------------------- line 152 -----------------------
        disp('Incorrect file format');
        disp('Incorrect file format')
----------------------- line 153 -----------------------
        return;
        return
----------------------- line 154 -----------------------
    end
    }
----------------------- line 155 -----------------------
    data = double(c.data); rows = c.rows; alleleCodes = c.alleleCodes;
    data <- double(c.data) rows <- c.rows alleleCodes <- c.alleleCodes
----------------------- line 156 -----------------------
    noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
    noalle <- c.noalle adjprior <- c.adjprior priorTerm <- c.priorTerm
----------------------- line 157 -----------------------
    dist = c.dist; Z = c.Z; popnames = c.popnames; rowsFromInd = c.rowsFromInd;
    dist <- c.dist Z <- c.Z popnames <- c.popnames rowsFromInd <- c.rowsFromInd
----------------------- line 158 -----------------------
    clear c;
    clear c
----------------------- line 159 -----------------------
end
}
----------------------- line 160 -----------------------


----------------------- line 161 -----------------------
c.data=data; c.rows = rows; c.alleleCodes = alleleCodes;
c.data=data c.rows <- rows c.alleleCodes <- alleleCodes
----------------------- line 162 -----------------------
c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
c.noalle <- noalle c.adjprior <- adjprior c.priorTerm <- priorTerm
----------------------- line 163 -----------------------
c.dist=dist; c.Z=Z; c.rowsFromInd = rowsFromInd;
c.dist=dist c.Z=Z c.rowsFromInd <- rowsFromInd
----------------------- line 164 -----------------------


----------------------- line 165 -----------------------
% partition compare
#  partition compare
----------------------- line 166 -----------------------
if ~isempty(partitionCompare)
if (!is.null(partitionCompare)) {
----------------------- line 167 -----------------------
    nsamplingunits = size(rows,1);
    nsamplingunits <- size(rows, 1)
----------------------- line 168 -----------------------
    partitions = partitionCompare.partitions;
    partitions <- partitionCompare.partitions
----------------------- line 169 -----------------------
    npartitions = size(partitions,2);
    npartitions <- size(partitions, 2)
----------------------- line 170 -----------------------
    partitionLogml = zeros(1,npartitions);
    partitionLogml <- zeros(1, npartitions)
----------------------- line 171 -----------------------
    for i = 1:npartitions
    for (i in 1:npartitions) {
----------------------- line 172 -----------------------
        % number of unique partition lables
       # number of unique partition lables
----------------------- line 173 -----------------------
        npops = length(unique(partitions(:,i)));
        npops <- length(unique(partitions[, i]))
----------------------- line 174 -----------------------
        try
        try
----------------------- line 175 -----------------------
            partitionInd = zeros(rows(end),1);
            partitionInd <- zeros(rows(end), 1)
----------------------- line 176 -----------------------
            partitionSample = partitions(:,i);
            partitionSample <- partitions[, i]
----------------------- line 177 -----------------------
            for j = 1: nsamplingunits
            for (j in 1:) { nsamplingunits
----------------------- line 178 -----------------------
                partitionInd([c.rows(j,1):c.rows(j,2)]) = partitionSample(j);
                partitionInd([c.rows(j, 1):c.rows[j, 2)]] <- partitionSample(j)
----------------------- line 179 -----------------------
            end
            }
----------------------- line 180 -----------------------
            partitionLogml(i) = ...
            partitionLogml[i] <-
----------------------- line 181 -----------------------
                initialCounts(partitionInd, data(:,1:end-1), npops, c.rows, noalle, adjprior);
                initialCounts(partitionInd, data[, 1:end - 1], npops, c.rows, noalle, adjprior)
----------------------- line 182 -----------------------
        catch
        catch
----------------------- line 183 -----------------------
           disp('*** ERROR: unmatched data.');
           disp(' * ** ERROR: unmatched data.')
----------------------- line 184 -----------------------
           return
           return
----------------------- line 185 -----------------------
        end
        }
----------------------- line 186 -----------------------
    end
    }
----------------------- line 187 -----------------------
    % return the logml result
   # return the logml result
----------------------- line 188 -----------------------
    partitionCompare.logmls = partitionLogml;
    partitionCompare.logmls <- partitionLogml
----------------------- line 189 -----------------------
    set(h1, 'userdata', partitionCompare);
    set(h1, 'userdata', partitionCompare)
----------------------- line 190 -----------------------
    return
    return
----------------------- line 191 -----------------------
end
}
----------------------- line 192 -----------------------


----------------------- line 193 -----------------------
if fixedK
if (fixedK) {
----------------------- line 194 -----------------------
    [logml, npops, partitionSummary]=indMix_fixK(c);
    [logml, npops, partitionSummary]=indMix_fixK(c)
----------------------- line 195 -----------------------
else
} else {
----------------------- line 196 -----------------------
    [logml, npops, partitionSummary]=indMix(c);
    [logml, npops, partitionSummary]=indMix(c)
----------------------- line 197 -----------------------
end
}
----------------------- line 198 -----------------------


----------------------- line 199 -----------------------
if logml==1
if (logml == 1) {
----------------------- line 200 -----------------------
    return
    return
----------------------- line 201 -----------------------
end
}
----------------------- line 202 -----------------------


----------------------- line 203 -----------------------
data = data(:,1:end-1);
data <- data[, 1:end - 1]
----------------------- line 204 -----------------------
viewPopMixPartition(PARTITION, rows, popnames);
viewPopMixPartition(PARTITION, rows, popnames)
----------------------- line 205 -----------------------
%npops = poistaTyhjatPopulaatiot(npops);
# npops <- poistaTyhjatPopulaatiot(npops)
----------------------- line 206 -----------------------
%POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
# POP_LOGML <- computePopulationLogml(1:npops, adjprior, priorTerm)
----------------------- line 207 -----------------------


----------------------- line 208 -----------------------
h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
h0 <- findobj('Tag', 'filename1_text')  inp <- get(h0, 'String')
----------------------- line 209 -----------------------
h0 = findobj('Tag','filename2_text');
h0 <- findobj('Tag', 'filename2_text')
----------------------- line 210 -----------------------
outp = get(h0,'String');
outp <- get(h0, 'String')
----------------------- line 211 -----------------------
changesInLogml = writeMixtureInfoPop(logml, rows, data, adjprior, priorTerm, ...
changesInLogml <- writeMixtureInfoPop(logml, rows, data, adjprior, priorTerm,
----------------------- line 212 -----------------------
    outp,inp,partitionSummary, popnames, fixedK);
    outp, inp, partitionSummary, popnames, fixedK)
----------------------- line 213 -----------------------


----------------------- line 214 -----------------------
talle = questdlg(['Do you want to save the mixture populations ' ...
talle <- questdlg(['Do you want to save the mixture populations '
----------------------- line 215 -----------------------
    'so that you can use them later in admixture analysis?'], ...
    'so that you can use them later in admixture analysis?'],
----------------------- line 216 -----------------------
    'Save results?','Yes','No','Yes');
    'Save results?', 'Yes', 'No', 'Yes')
----------------------- line 217 -----------------------
if isequal(talle,'Yes')
if (isequal(talle, 'Yes')) {
----------------------- line 218 -----------------------
    waitALittle;
    waitALittle
----------------------- line 219 -----------------------
    [filename, pathname] = uiputfile('*.mat','Save results as');
    [filename, pathname] = uiputfile(' * .mat', 'Save results as')
----------------------- line 220 -----------------------


----------------------- line 221 -----------------------
    if (filename == 0) & (pathname == 0)
    if ((filename == 0) & (pathname == 0)) {
----------------------- line 222 -----------------------
        % Cancel was pressed
       # Cancel was pressed
----------------------- line 223 -----------------------
        return
        return
----------------------- line 224 -----------------------
    else % copy 'baps4_output.baps' into the text file with the same name.
    else# copy 'baps4_output.baps' into the text file with the same name.
----------------------- line 225 -----------------------
        if exist('baps4_output.baps','file')
        if (exist('baps4_output.baps', 'file')) {
----------------------- line 226 -----------------------
            copyfile('baps4_output.baps',[pathname filename '.txt'])
            copyfile('baps4_output.baps', [pathname filename '.txt'])
----------------------- line 227 -----------------------
            delete('baps4_output.baps')
            delete('baps4_output.baps')
----------------------- line 228 -----------------------
        end
        }
----------------------- line 229 -----------------------
    end
    }
----------------------- line 230 -----------------------


----------------------- line 231 -----------------------
    if rowsFromInd==0
    if (rowsFromInd == 0) {
----------------------- line 232 -----------------------
        %K�ytettiin BAPS-formaattia, eik?rowsFromInd ole tunnettu.
       # K < e4 > ytettiin BAPS - formaattia, eik?rowsFromInd ole tunnettu.
----------------------- line 233 -----------------------
        [popnames, rowsFromInd] = findOutRowsFromInd(popnames, rows);
        [popnames, rowsFromInd] = findOutRowsFromInd(popnames, rows)
----------------------- line 234 -----------------------
    end
    }
----------------------- line 235 -----------------------


----------------------- line 236 -----------------------
    groupPartition = PARTITION;
    groupPartition <- PARTITION
----------------------- line 237 -----------------------


----------------------- line 238 -----------------------
    fiksaaPartitioYksiloTasolle(rows, rowsFromInd);
    fiksaaPartitioYksiloTasolle(rows, rowsFromInd)
----------------------- line 239 -----------------------


----------------------- line 240 -----------------------
    c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
    c.PARTITION <- PARTITION c.COUNTS <- COUNTS c.SUMCOUNTS <- SUMCOUNTS
----------------------- line 241 -----------------------
    c.alleleCodes = alleleCodes; c.adjprior = adjprior;
    c.alleleCodes <- alleleCodes c.adjprior <- adjprior
----------------------- line 242 -----------------------
    c.rowsFromInd = rowsFromInd; c.popnames = popnames;
    c.rowsFromInd <- rowsFromInd c.popnames <- popnames
----------------------- line 243 -----------------------
    c.data = data; c.npops = npops; c.noalle = noalle;
    c.data <- data c.npops <- npops c.noalle <- noalle
----------------------- line 244 -----------------------
    c.mixtureType = 'popMix'; c.groupPartition = groupPartition;
    c.mixtureType = 'popMix' c.groupPartition <- groupPartition
----------------------- line 245 -----------------------
    c.rows = rows; c.logml = logml; c.changesInLogml = changesInLogml;
    c.rows <- rows c.logml <- logml c.changesInLogml <- changesInLogml
----------------------- line 246 -----------------------
%     save([pathname filename], 'c');
#      save([pathname filename], 'c')
----------------------- line 247 -----------------------
    save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
    save([pathname filename], 'c', ' - v7.3')# added by Lu Cheng, 08.06.2012
----------------------- line 248 -----------------------
else
} else {
----------------------- line 249 -----------------------
    if exist('baps4_output.baps','file')
    if (exist('baps4_output.baps', 'file')) {
----------------------- line 250 -----------------------
        delete('baps4_output.baps')
        delete('baps4_output.baps')
----------------------- line 251 -----------------------
    end
    }
----------------------- line 252 -----------------------
end
}
----------------------- line 253 -----------------------
NA
        return(function greedyPopMix)
}
