function trainedMix

% LASKENNAN ALKUARVOJEN Mï¿½ï¿½RITTï¿½MINEN

global SCRIPT_MODE;

if isempty(SCRIPT_MODE)
    SCRIPT_MODE = false;
end

if SCRIPT_MODE
    input_type = 'MLST-format';
else
    input_type = questdlg('Specify the format of your data: ',...
                          'Specify Data Format', ...
                          'MLST-format', 'GenePop-format','MLST-format');
end

switch input_type
    case 'MLST-format'
        disp('MLST-format');
        processMLST        
    case 'GenePop-format'
        disp('GenePop-format');
        processGenePop
end

%--------------------------------------------------------------------------

function processMLST
% note that this version only works for windows with Excel installed
% Lu Cheng, 02.02.2010
% lu.cheng@helsinki.fi

% added by Lu Cheng, 11.03.2010
global SCRIPT_MODE;
global PARAMETERS;
if isempty(SCRIPT_MODE)
    SCRIPT_MODE = false;
end
% ----------

tmp_train_file = 'tmp8972_train.xls';
if exist(tmp_train_file,'file')==2
    delete(tmp_train_file);
end

%% process both the training data and test data

% Format of the training excel file
% column 1: sample ID
% column 2: cluster label of each sample, an integer from 1 to K
% column 3-n: sequences of each gene
format1 = 'MS EXCEL FORMAT';
format2 = 'PREPROCESSED FORMAT';

if SCRIPT_MODE
    if isequal(PARAMETERS.train_file_format,'.xls')
        input_type = format1;
    elseif isequal(PARAMETERS.train_file_format,'.mat')
        input_type = format2;
    end
else
    input_type = questdlg('Specify the format of your training data: ',...
    'Specify Data Format', format1, format2, format1);
end

switch input_type
    case format1
        
        if SCRIPT_MODE
            trained_file = PARAMETERS.train_file_name;
        else
            [filename, pathname] = uigetfile('*.xls', strcat('Load training data in',' ',format1));
            if filename==0
                return;
            end
            trained_file = strcat(pathname,filename);
        end
                
        [A B] = xlsread(trained_file);
        if size(B,1) == length(A)+1
            B(2:end,1) = num2cell(A(:,1));
        else
            B(:,1) = num2cell(A(:,1));
        end

        train_xls = B(:,[1 3:end]);
        cluster_labels = A(:,2);
        
        % the unique labels should be tightly from 1 to K
        % added by Lu Cheng, 22.06.2010
        unique_labels = unique(cluster_labels);
        if max(unique_labels)~=length(unique_labels)
            error('The cluster labels are wrong, should be from 1 to %s !', num2str(length(unique_labels)));
        end
        
        xlswrite(tmp_train_file,train_xls);
        clear A B trained_file unique_labels

        c_train = preprocessXLS(tmp_train_file);
        c_train.cluster_labels = cluster_labels;
        delete(tmp_train_file);
        
        if SCRIPT_MODE
            save_preproc = PARAMETERS.save_prepro_train_data;
        else
            save_preproc = questdlg('Do you wish to save the pre-processed training data?',...
                                        'Save pre-processed data?',...
                                   'Yes','No','Yes');
        end
        
        if isequal(save_preproc,'Yes');
            if SCRIPT_MODE
%                save(PARAMETERS.train_prepro_file,'c_train');
               save(PARAMETERS.train_prepro_file,'c_train', '-v7.3'); % added by Lu Cheng, 08.06.2012
            else
               [filename, pathname] = uiputfile('*.mat','Save pre-processed training data as');
               if (sum(filename)==0) || (sum(pathname)==0)
                   % do nothing
               else
%                    save(strcat(pathname,filename,'.mat'),'c_train');
                   save(strcat(pathname,filename,'.mat'),'c_train','-v7.3'); % added by Lu Cheng, 08.06.2012
               end
           end
        end;
        
    case format2
        disp(format2);
        
        if SCRIPT_MODE
            trained_file = PARAMETERS.train_file_name;
        else
            [filename, pathname] = uigetfile('*.mat', strcat('Load training data in',' ',format2));
            if filename==0
                return;
            end
            trained_file = strcat(pathname,filename);
        end
        
        clear c_train
        load('-mat',trained_file);
    
    otherwise
        return;
end

%% process with test data

if SCRIPT_MODE
    if isequal(PARAMETERS.test_file_format,'.xls')
        input_type = format1;
    elseif isequal(PARAMETERS.test_file_format,'.mat')
        input_type = format2;
    end
else
    input_type = questdlg('Specify the format of your test data: ',...
    'Specify Data Format', format1, format2, format1);
end

switch input_type
    case format1
        if SCRIPT_MODE
            test_file = PARAMETERS.test_file_name;
        else
            [filename, pathname] = uigetfile('*.xls', 'Load test data (unlabeled) in MLST-format');
            if filename==0
                return;
            end
            test_file = strcat(pathname,filename);
        end
        c_test = preprocessXLS(test_file,c_train);
        
        if SCRIPT_MODE
            save_preproc = PARAMETERS.save_prepro_test_data;
        else
            save_preproc = questdlg('Do you wish to save the pre-processed test data?',...
                                        'Save pre-processed data?','Yes','No','Yes');
        end
        
        if isequal(save_preproc,'Yes');
            if SCRIPT_MODE
%                 save(PARAMETERS.test_prepro_file,'c_test');
                save(PARAMETERS.test_prepro_file,'c_test','-v7.3'); % added by Lu Cheng, 08.06.2012
            else
               [filename, pathname] = uiputfile('*.mat','Save pre-processed test data as');
               if (sum(filename)==0) || (sum(pathname)==0)
                   % do nothing
               else
%                    save(strcat(pathname,filename,'.mat'),'c_test');
                   save(strcat(pathname,filename,'.mat'),'c_test','-v7.3'); % added by Lu Cheng, 08.06.2012
               end
            end
        end;

    case format2
        if SCRIPT_MODE
            test_file = PARAMETERS.test_file_name;
        else
            [filename, pathname] = uigetfile('*.mat', cat(2,'Load test data (unlabeled) in ',format2));
            if filename==0
                return;
            end
            test_file = strcat(pathname,filename);
        end
        load('-mat',test_file,'c_test');
        
    otherwise
        return;
end

%% compare the preprocessed training and test data and further steps

semi_linkageMixture_speed(c_train, c_test);

%--------------------------------------------------------------------------
function processGenePop
global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global ADJPRIOR;
global PRIORTERM;
global SUMPRIOR;
global LOGDIFF;
clearGlobalVars;

[filename, pathname] = uigetfile('*.txt', 'Load prior data in GenePop-format');
if filename==0
    return;
end
kunnossa = testaaGenePopData([pathname filename]);
if kunnossa==0
    return
end

waitALittle;
[filename2, pathname2] = uigetfile('*.txt', 'Load sampling units in GenePop-format');
if filename2==0
    return;
end
kunnossa = testaaGenePopData([pathname2 filename2]);
if kunnossa==0
    return
end
clear kunnossa;

[pData, pNames, pIndNames]=lueGenePopDataPop([pathname filename]);
[suData, suNames, suIndNames] = lueGenePopDataPop([pathname2 filename2]);
if size(pData,2) ~= size(suData,2)
    disp('Incorrect input');
    return
end
inp = [filename ' & ' filename2];
h0 = findobj('Tag','filename1_text');
set(h0,'String',inp);
clear h0; clear inp;
clear filename; clear filename2; clear pathname; clear pathname2;

[alleleCodes, noalle, suData, pData] = examineAlleles(suData, pData);

rows = initializeRows(suData);   % Samplin unit:ien rivit kertova muuttuja.
rowsFromInd = 2;  %Tiedetï¿½ï¿½n GenePop:in tapauksessa.
data = suData(:,1:end-1);   %Klusteroitavat "yksilï¿½t"
priorLastCol = pData(:,end); 
priorPartition = priorLastCol(1:rowsFromInd:end);  % Prioriyksilï¿½iden partitio
clear suData; clear priorLastCol; %Ei tarvita. Kai...?

npopstext = [];
ready = false;
teksti = 'Input upper bound to the number of populations (possibly multiple values): ';
while ready == false
    npopstextExtra = inputdlg(teksti ,...
        'Input maximum number of populations',1,{'20'});
    if isempty(npopstextExtra)  % Painettu Cancel:ia
        return
    end
    npopstextExtra = npopstextExtra{1};
    if length(npopstextExtra)>=255
        npopstextExtra = npopstextExtra(1:255);
        npopstext = [npopstext ' ' npopstextExtra];
        teksti = 'The input field length limit (255 characters) was reached. Input more values: ';
    else
        npopstext = [npopstext ' ' npopstextExtra];
        ready = true;
    end
end
clear ready; clear teksti;
if isempty(npopstext) | length(npopstext)==1
    return
else
    npopsTaulu = str2num(npopstext);
    clear npopstext;
    if length(npopsTaulu)<1
        disp('Incorrect input');
        return
    end
    if any(npopsTaulu < size(pNames,1))
        disp('Incorrect input');
        return
    end
end

nruns = length(npopsTaulu);

logmlBest = -1e50;
partitionSummary = -1e50*ones(30,2);  % Tiedot 30 parhaasta partitiosta (npops ja logml)
partitionSummary(:,1) = zeros(30,1);
worstLogml = -1e50; worstIndex = 1;
Z = [];
for run = 1:nruns 
    npops = npopsTaulu(run);
    dispLine;
    disp(['Run ' num2str(run) '/' num2str(nruns) ...
        ', maximum number of populations ' num2str(npops) '.']); 
    
    disp(['Simulation started with ' num2str(npops) ' initial populations.']);
    adjprior = computePriors(pData, npops, noalle);  %adjprior on yhden populaation, jossa ei havaintoja.
    COUNTS = zeros(size(ADJPRIOR));
    SUMCOUNTS = zeros(size(SUMPRIOR));
    POP_LOGML = zeros(npops,1);

    POP_LOGML = computePopulationLogml(1:npops);
    logml = initialPopCounts(data, npops, rows, noalle);  %Alustetaan COUNTS, PARTITION ...
    
    if isempty(Z) % Lasketaan vain ensimmï¿½isellï¿?kierroksella.
        if size(rows,1)==1
            Z = [];
            dist = [];
        else
            [Z,dist] = getPopDistancesByKL(data, rows, noalle, adjprior);  %Lasketaan sampling unit:ien vï¿½liset etï¿½isyydet.    
        end
    end
        
    if logml>worstLogml
        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
        if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
    end
    
    % PARHAAN MIXTURE-PARTITION ETSIMINEN
    nRoundTypes = 7;
    kokeiltu = zeros(nRoundTypes, 1);

    roundTypes = [1 1];  %Ykkï¿½svaiheen sykli kahteen kertaan.   
    ready = 0; vaihe = 1;
    ninds = length(PARTITION);  % num of sampling units
    LOGDIFF = repmat(-Inf,ninds,npops);

    disp(' ');
    while ready ~= 1
	    muutoksia = 0;
	
        disp(['Performing steps: ' num2str(roundTypes)]);
    
        for n = 1:length(roundTypes)
        
            round = roundTypes(n);
            kivaluku=0;
            
            if kokeiltu(round) == 1
        
            elseif round==0 | round==1   %Yksilï¿½n siirtï¿½minen toiseen populaatioon.
                inds = 1:ninds;
                aputaulu = [inds' rand(ninds,1)];
                aputaulu = sortrows(aputaulu,2);
                inds = aputaulu(:,1)';
            
                muutosNyt = 0;
                for ind = inds
                    i1 = PARTITION(ind);
                    [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
                        data);
                
                    if round==1, [maxMuutos, i2] = max(muutokset);
                    end
                
                    if (i1~=i2 & maxMuutos>1e-5)
                        % Tapahtui muutos
                        if muutosNyt == 0
                            disp('Action 1');
                            muutosNyt = 1;
                            kokeiltu = zeros(nRoundTypes,1);
                        end
                        muutoksia = 1;
                        kivaluku = kivaluku+1;
                        updateGlobalVariables(ind, i2, diffInCounts);
                        logml = logml+maxMuutos;
                        
                        if logml>worstLogml
                            [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                            if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                        end    
                    end
                end
                if muutosNyt == 0
                   kokeiltu(round) = 1;
                end 
            
            elseif round==2 & ~isempty(dist)   %Populaation yhdistï¿½minen toiseen.
                maxMuutos = 0;
                for pop = 1:npops
                    [muutokset, diffInCounts] = laskeMuutokset2(pop, rows, ...
                        data);
                    [isoin, indeksi] = max(muutokset);
                    if isoin>maxMuutos
                        maxMuutos = isoin;
                        i1 = pop;
                        i2 = indeksi;
                        diffInCountsBest = diffInCounts;
                    end
                end
            
                if maxMuutos>1e-5
                    muutoksia = 1;
                    disp('Action 2');
                    kokeiltu = zeros(nRoundTypes,1);
                    updateGlobalVariables2(i1,i2, diffInCountsBest);
                    logml = logml + maxMuutos;
                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                        if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                    end
                else
                    kokeiltu(round) = 1;
                end
                        
            
            elseif (round==3 | round==4) & ~isempty(dist)%Populaation jakaminen osiin.
                maxMuutos = 0;
                ninds = size(rows,1);
                for pop = 1:npops
                    inds2 = find(PARTITION==pop);
                    ninds2 = length(inds2);
                    if ninds2>2
                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = linkage(dist2');
                        if round==3
                            npops2 = min(20, floor(ninds2 / 5));
                        elseif round==4
                            npops2 = 2;  %Moneenko osaan jaetaan
                        end
                        T2 = cluster_own(Z2, npops2);
                        muutokset = laskeMuutokset3(T2, inds2, rows, data, pop);
                        [isoin, indeksi] = max(muutokset(1:end));
                        if isoin>maxMuutos
                            maxMuutos = isoin;
                            muuttuvaPop2 = rem(indeksi,npops2);
                            if muuttuvaPop2==0, muuttuvaPop2 = npops2; end
                            muuttuvat = inds2(find(T2==muuttuvaPop2));
                            i2 = ceil(indeksi/npops2);
                        end
                    end
                end
                if maxMuutos>1e-5
                    muutoksia = 1;
                    disp(['Action ' num2str(round)]);
                    kokeiltu = zeros(nRoundTypes,1);
               	    %rows = computeRows(rowsFromInd, muuttuvat, length(muuttuvat));
                    rivit = [];
                    for ind = muuttuvat
                        lisa = rows(ind,1):rows(ind,2);
                        rivit = [rivit; lisa'];
                        %rivit = [rivit; rows(ind)'];
                    end
                    diffInCounts = computeDiffInCounts(rivit', size(COUNTS,1), ...
                    size(COUNTS,2), data);
                    i1 = PARTITION(muuttuvat(1));
                    updateGlobalVariables3(muuttuvat, diffInCounts, i2);
                    logml = logml + maxMuutos;
                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                        if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                    end
                else
                    kokeiltu(round)=1;
                end
            
            
            elseif round == 5 & ~isempty(dist)
                % Kï¿½y lï¿½pi populaatioita.
                % Yritï¿?poistaa niistï¿?yksilï¿½itï¿?yksi
                % kerrallaan. Lopeta heti, kun jonkin 
                % yksilï¿½iden joukon poistaminen jostain
                % populaatiosta aiheuttaa positiivisen 
                % muutoksen logml:ï¿½ï¿½n.
            
                pop=0;
                muutettu = 0;
                poplogml = POP_LOGML;    partition = PARTITION;
                counts = COUNTS;         sumcounts = SUMCOUNTS;
                logdiff = LOGDIFF;
            
                while (pop < npops & muutettu == 0)
                    pop = pop+1;
                    totalMuutos = 0;
                    inds = find(PARTITION==pop)';
                    inds = returnInOrder(inds, pop, rows, data);
               
                    i=0;
                    while (length(inds)>0 & i<length(inds) & totalMuutos<1e-5) % Lopetetaankun totalMuutos > 0
                        i = i+1;
                        ind = inds(i);
                        [muutokset, diffInCounts] = laskeMuutokset(ind, rows, data);
                        muutokset(pop) = -1e50;   % Varmasti ei suurin!!!
                        [maxMuutos, i2] = max(muutokset);
                        updateGlobalVariables(ind, i2, diffInCounts);
                        totalMuutos = totalMuutos+maxMuutos;
                        logml = logml+maxMuutos;
                    end
               
                    if totalMuutos>1e-5
                        disp('action 5');
                        muutettu=1;
                        kokeiltu = zeros(nRoundTypes,1);
                        muutoksia = 1;  % Ulompi kirjanpito.
                        if logml>worstLogml
                            [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                            if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                        end    
                    else
                        % Missï¿½ï¿½n vaiheessa tila ei parantunut.
                        % Perutaan kaikki muutokset.
                        PARTITION = partition;
                        SUMCOUNTS = sumcounts;
                        POP_LOGML = poplogml;
                        COUNTS = counts;
                        LOGDIFF = logdiff;
                        logml = logml - totalMuutos;
                        kokeiltu(round)=1;
                    end
                end
                clear partition; clear sumcounts; clear counts; clear poplogml;
            end
        end
    
        if muutoksia == 0
            if vaihe==1
                vaihe = 2;
            elseif vaihe==2
                vaihe = 3;
            elseif vaihe==3
                ready = 1;
            end
        else
            muutoksia = 0;
        end
    
        if ready==0
            if vaihe==1
                roundTypes=[1];
            elseif vaihe==2
                roundTypes=[2 1];
            elseif vaihe==3
                roundTypes=[5 4 3 1 2];
            end
        end         
    end


    % TALLENNETAAN

    prioriPopLkm = size(pNames,1);
    npops = poistaTyhjatPopulaatiot(prioriPopLkm);
    POP_LOGML = computePopulationLogml(1:npops);

    n_clust_with_su = length(unique(PARTITION));
    disp(['Found partition with sampling units in ' num2str(n_clust_with_su) ' clusters.']);
    disp(['Log(ml) = ' num2str(logml)]);
    disp(' ');
    
    if logml>logmlBest
        % Pï¿½ivitetï¿½ï¿½n parasta lï¿½ydettyï¿?partitiota.
        logmlBest = logml;
        npopsBest = npops;
        partitionBest = PARTITION;
        countsBest = COUNTS;
        sumCountsBest = SUMCOUNTS;
        pop_logmlBest = POP_LOGML;
        adjPriorBest = ADJPRIOR;
        priorTermBest = PRIORTERM;
        sumPriorBest = SUMPRIOR;
        logdiffbest = LOGDIFF;
    end    
end

logml = logmlBest;
npops = npopsBest;
PARTITION = partitionBest;
COUNTS = countsBest;
SUMCOUNTS = sumCountsBest;
POP_LOGML = pop_logmlBest;
ADJPRIOR = adjPriorBest;
PRIORTERM = priorTermBest;
SUMPRIOR = sumPriorBest;
LOGDIFF = logdiffbest;

h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
h0 = findobj('Tag','filename2_text');  outp = get(h0,'String');
writeTrainedMixtureInfo(logml, rows, data, outp, inp, ...
    suIndNames, suNames, pIndNames, pNames, partitionSummary);

fiksaaPartitioYksiloTasolle(rows, rowsFromInd);
[data, popnames] = muokkaaMuuttujat(adjprior, rowsFromInd, ...
    pNames, suNames, priorPartition, pData, data);
viewMixPartition(PARTITION, popnames);

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

    c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
    c.alleleCodes = alleleCodes; c.adjprior = adjprior;
    c.rowsFromInd = rowsFromInd; c.popnames = popnames;
    c.data = data; c.npops = npops; c.noalle = noalle;
    c.mixtureType = 'trained';
%     save([pathname filename], 'c');
    save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
else
    if exist('baps4_output.baps','file')
        delete('baps4_output.baps')
    end
end

%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedetï¿½ï¿½n, ettï¿?annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ssï¿?ei vielï¿?ole
% annettua logml arvoa, niin lisï¿½tï¿½ï¿½n worstIndex:in kohtaan uusi logml ja
% nykyistï¿?partitiota vastaava nclusters:in arvo. Muutoin ei tehdï¿?mitï¿½ï¿½n.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lï¿½ydetty partitio ei ole vielï¿?kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end


%--------------------------------------------------------------------------



function [data, popnames] = muokkaaMuuttujat(adjprior, rowsFromInd, ...
    pNames, suNames, priorPartition, pData, data)
% Muokkaa kaikki tarvittavat muuttujat mixture result-file
% muotoisiksi.

global PARTITION; global COUNTS;   
global SUMCOUNTS; global ADJPRIOR;
nloci = size(data,2);
npops = size(COUNTS, 3);

data = [pData(:,1:nloci) ; data];
PARTITION = [priorPartition; PARTITION];
priorCounts = ADJPRIOR-repmat(adjprior, [1 1 npops]);
COUNTS = COUNTS+priorCounts;
SUMCOUNTS = (squeeze(sum(COUNTS)))';

priorNinds = length(priorPartition);
for k = 1:size(suNames,1)
    suNames{k,2} = suNames{k,2} + priorNinds;
end
popnames = [pNames; suNames];


%-------------------------------------------------------------------------
    
    
function [alleleCodes, noalle, suData, pData] = examineAlleles(suData, pData)
% Poistetaan nollat molemmista datoista. Selvitetï¿½ï¿½n noalle ja
% alleleCodes ja muutetaan molemmat datat vastaamaan alleleCodes:ia.
% Tï¿½ssï¿?vaiheessa datojen viimeinen sarake kertoo yksikï¿½n, jolle
% rivi kuuluu.
data = [pData; suData];
nrows_prior = size(pData,1);
nloci = size(suData,2)-1;

dataApu = data(:,1:nloci); %poistetaan nollat
nollat = find(dataApu==0); 
if ~isempty(nollat)
   isoinAlleeli = max(max(dataApu));
   dataApu(nollat) = isoinAlleeli+1;
   data(:,1:nloci) = dataApu;
end
dataApu = []; nollat = []; isoinAlleeli = [];

noalle=zeros(1,nloci); %selvitetï¿½ï¿½n noalle
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
end

alleleCodes = zeros(max(noalle),nloci); %selvitetï¿½ï¿½n alleleCodes
for i=1:nloci  
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
end

for loc = 1:nloci   %muutetaan alleelien koodit vastaamaan alleleCodes:ia
    for all = 1:noalle(loc)
        data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;
pData = data(1:nrows_prior , :);
suData = data(nrows_prior+1:end , :);


%----------------------------------------------------------------------

function adjprior = computePriors(pData, npops, noalle)
global ADJPRIOR;
global SUMPRIOR;
global PRIORTERM;
nloci = size(pData,2)-1;
max_noalle = max(noalle);
ADJPRIOR = zeros(max_noalle, nloci, npops);
PRIORTERM = zeros(npops, 1);
SUMPRIOR = zeros(npops, nloci);

adjprior = zeros(max_noalle,nloci);
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
end

data = pData(:, 1:nloci);
for i = 1:npops
    rivit = find(pData(:,end) == i)'; %Pitï¿½ï¿½ olla vaakavektori.
    if ~isempty(rivit)
        diffInCounts = computeDiffInCounts(rivit, max_noalle, nloci, data);
        ADJPRIOR(:,:,i) = diffInCounts;
    end
    ADJPRIOR(:,:,i) = ADJPRIOR(:,:,i) + adjprior;
    for j=1:nloci
        SUMPRIOR(i,j) = sum(squeeze(ADJPRIOR(1:noalle(j), j , i)));
        PRIORTERM(i) = PRIORTERM(i)+gammaln(SUMPRIOR(i,j));
        PRIORTERM(i) = PRIORTERM(i)-sum(gammaln(squeeze(ADJPRIOR(1:noalle(j),j,i))));
    end
end


%--------------------------------------------------------------

function rows = initializeRows(data)
% Lasketaan rows-muuttuja. Tï¿½ssï¿?vaiheessa datan
% viimeisessï¿?sarakkeessa on vielï¿?yksikï¿½n kertova
% indeksi.
nind = max(data(:,end));
rows = zeros(nind,2);
for i=1:nind
    rivit = find(data(:,end)==i)';
    rows(i,1) = min(rivit);
    rows(i,2) = max(rivit);
end


%----------------------------------------------------------------


function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];
global ADJPRIOR; ADJPRIOR = [];
global PRIORTERM; PRIORTERM = [];
global SUMPRIOR; SUMPRIOR = [];
global LOGDIFF; LOGDIFF = [];


%--------------------------------------------------------------------

function [Z,distances] = getPopDistancesByKL(data, rows, noalle, adjprior)
% Laskee populaatioille etï¿½isyydet
% kï¿½yttï¿½en KL-divergenssiï¿?

npops = size(rows,1);   %Samplin unit:tien lkm
nloci=size(data,2);
maxnoalle = max(noalle);
counts = zeros(maxnoalle,nloci,npops);   % Tilapï¿½istï¿?kï¿½yttï¿½ï¿½ varten
sumcounts = zeros(npops,nloci);

for i=1:npops
    for j=1:nloci
        i_rivit = rows(i,1):rows(i,2);
        havainnotLokuksessa = find(data(i_rivit,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(i_rivit,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

distances = zeros(nchoosek(npops,2),1);

d = zeros(maxnoalle, nloci, npops);
prior = adjprior;
prior(find(prior==1))=0;
nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yhtï¿?alleelia.
prior(1,nollia)=1;

for pop1 = 1:npops
    d(:,:,pop1) = (squeeze(counts(:,:,pop1))+prior) ./ repmat(sum(squeeze(counts(:,:,pop1))+prior),maxnoalle,1);
    %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
end
pointer = 1;
for pop1 = 1:npops-1
    for pop2 = pop1+1:npops
        dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
        div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
        div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
        div = (div12+div21)/2;
        distances(pointer) = div;
        pointer = pointer+1;    
    end
end
Z=linkage(distances');

%--------------------------------------------------------------------------


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


%-----------------------------------------------------------------------


function logml = initialPopCounts(data, npops, rows, noalle)

global COUNTS;
global SUMCOUNTS;
global PARTITION;
global POP_LOGML;
global ADJPRIOR;
global SUMPRIOR;
nloci=size(data,2);
ninds = size(rows,1);
COUNTS = zeros(max(noalle),nloci,npops);
SUMCOUNTS = zeros(npops,nloci);
PARTITION = zeros(1,ninds);

inds = 1:ninds;
aputaulu = [inds' rand(ninds,1)];
aputaulu = sortrows(aputaulu,2);
inds = aputaulu(:,1)';
%omaPartitio = 1:6;    %POIS!!!!!!!!
%omaPartitio = omaPartitio';
%omaPartitio = omaPartitio(:,ones(30,1));
%omaPartitio = omaPartitio';
%omaPartitio = omaPartitio(:);   %POIS
%keyboard;
for ind = inds  % Sijoitetaan yksilï¿½t yksi kerrallaan.
    [muutokset, diffInCounts] = ...
        laskePrioriMuutokset(ind, rows, data);
    [maxMuutos, i2] = max(muutokset);
    %i2 = omaPartitio(ind)    %POIS
    PARTITION(ind) = i2;
    COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
    SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);
    POP_LOGML(i2) = computePopulationLogml(i2);
    
end

logml = laskeLoggis(COUNTS, SUMCOUNTS, ADJPRIOR, SUMPRIOR);

%keyboard;

%-----------------------------------------------------------------------


function loggis = laskeLoggis(counts, sumcounts, adjprior, sumprior)
npops = size(counts,3);

logml2 = sum(sum(sum(gammaln(counts+adjprior)))) ...
    - sum(sum(sum(gammaln(adjprior)))) ...
    - sum(sum(gammaln(sumcounts+sumprior))) ...
    + sum(sum(gammaln(sumprior)));
loggis = logml2;


%--------------------------------------------------------------------


function kunnossa = testaaGenePopData(tiedostonNimi)
% kunnossa == 0, jos data ei ole kelvollinen genePop data.
% Muussa tapauksessa kunnossa == 1.

kunnossa = 0;
fid = fopen(tiedostonNimi);
line1 = fgetl(fid);  %ensimmï¿½inen rivi
line2 = fgetl(fid);  %toinen rivi
line3 = fgetl(fid);  %kolmas

if (isequal(line1,-1) | isequal(line2,-1) | isequal(line3,-1))
    disp('Incorrect file format'); fclose(fid);
    return 
end
if (testaaPop(line1)==1 | testaaPop(line2)==1)
    disp('Incorrect file format'); fclose(fid);
    return
end
if testaaPop(line3)==1
    %2 rivi tï¿½llï¿½in lokusrivi
    nloci = rivinSisaltamienMjonojenLkm(line2);
    line4 = fgetl(fid);
    if isequal(line4,-1)
        disp('Incorrect file format'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin neljï¿?tï¿½ytyy sisï¿½ltï¿½ï¿½ pilkku.
        disp('Incorrect file format'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetï¿½ï¿½n, ettï¿?pysï¿½htyy
        pointer = pointer+1;
    end
    line4 = line4(pointer+1:end);  %pilkun jï¿½lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format'); fclose(fid);
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
        disp('Incorrect file format'); fclose(fid);
        return
    end
    nloci = lineNumb-2;
    line4 = fgetl(fid);  %Eka rivi pop sanan jï¿½lkeen
    if isequal(line4,-1)
        disp('Incorrect file format'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin tï¿½ytyy sisï¿½ltï¿½ï¿½ pilkku.
        disp('Incorrect file format'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetï¿½ï¿½n, ettï¿?pysï¿½htyy.
        pointer = pointer+1;
    end
 
    line4 = line4(pointer+1:end);  %pilkun jï¿½lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format'); fclose(fid);
        return
    end
end
kunnossa = 1;
fclose(fid);

%--------------------------------------------------------------------


function [data, popnames, indnames] = lueGenePopDataPop(tiedostonNimi)
% Data annetaan muodossa, jossa viimeinen sarake kertoo ryhmï¿½n.
% popnames on kuten ennenkin.

fid = fopen(tiedostonNimi);
line = fgetl(fid);  %ensimmï¿½inen rivi
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
indnames = cell(100,1);
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
        popnames{nimienLkm, 1} = {nimi};   %Nï¿½in se on greedyMix:issï¿½kin?!?
        popnames{nimienLkm, 2} = ninds;
        poimiNimi=0;
        
        data = addAlleles(data, ninds, line, divider);
        
        if ninds>size(indnames,1)
            indnames = [indnames; cell(100,1)];
        end
        indnames{ninds} = {nimi};
        
    elseif testaaPop(line)
        poimiNimi = 1;
        
    elseif line ~= -1
        ninds = ninds+1;
        nimi = lueNimi(line);
        data = addAlleles(data, ninds, line, divider);
        
        if ninds>size(indnames,1)
            indnames = [indnames; cell(100,1)];
        end
        indnames{ninds} = {nimi};
    end
end

indnames = indnames(1:ninds);
fclose(fid);
data = data(1:ninds*2,:);
popnames = popnames(1:nimienLkm,:);
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
% line on ensimmï¿½inen pop-sanan jï¿½lkeinen rivi
% Genepop-formaatissa olevasta datasta. funktio selvittï¿½ï¿½
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
% Palauttaa line:n sisï¿½ltï¿½mien mjonojen lukumï¿½ï¿½rï¿½n.
% Mjonojen vï¿½lissï¿?tï¿½ytyy olla vï¿½lilyï¿½nti.
count = 0;
pit = length(line);
tila = 0;    %0, jos odotetaan vï¿½lilyï¿½ntejï¿? 1 jos odotetaan muita merkkejï¿?
for i=1:pit
    merkki = line(i);
    if (isspace(merkki) & tila==0) 
        %Ei tehdï¿?mitï¿½ï¿½n.
    elseif (isspace(merkki) & tila==1)
        tila = 0;
    elseif (~isspace(merkki) & tila==0)
        tila = 1;
        count = count+1;
    elseif (~isspace(merkki) & tila==1)
        %Ei tehdï¿?mitï¿½ï¿½n
    end
end

%-------------------------------------------------------

function pal = testaaPop(rivi)
% pal=1, mikï¿½li rivi alkaa jollain seuraavista
% kirjainyhdistelmistï¿? Pop, pop, POP. Kaikissa muissa
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
% yksilï¿½ï¿½ ind vastaavat rivit. Yksilï¿½n alleelit
% luetaan genepop-formaatissa olevasta rivistï¿?
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

%------------------------------------------------------------------------------------


function popLogml = computePopulationLogml(pops)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on mï¿½ï¿½ritelty pops-muuttujalla.

global ADJPRIOR;
global PRIORTERM;
global SUMPRIOR;
global COUNTS;
global SUMCOUNTS;
x = size(COUNTS,1);
y = size(COUNTS,2);
z = length(pops);

%popLogml = ...
%    squeeze(sum(sum(reshape(...
%    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
%    ,[x y z])))) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;

popLogml = ...
    squeeze(sum(sum(reshape( ...
    gammaln(ADJPRIOR(:,:,pops) + COUNTS(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(SUMPRIOR(pops,:)+SUMCOUNTS(pops,:)),2) + PRIORTERM(pops);

%--------------------------------------------------------------------------


function [muutokset, diffInCounts] = ...
    laskePrioriMuutokset(ind, globalRows, data)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikï¿?olisi
% muutos logml:ssï¿? mikï¿½li yksilï¿?ind LISï¿½Tï¿½ï¿½N koriin i.

global COUNTS;  global SUMCOUNTS;
global POP_LOGML;
npops = size(COUNTS,3);
muutokset = zeros(npops,1);

rows = globalRows(ind,1):globalRows(ind,2);
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

i2 = [1:npops];

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops 1]);
new_i2_logml = computePopulationLogml(i2);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops 1]);

muutokset(i2) = new_i2_logml - POP_LOGML;


%--------------------------------------------------------------------------


function inds = returnInOrder(inds, pop, globalRows, data)
% Palauttaa yksilï¿½t sellaisessa jï¿½rjestyksessï¿? ettï¿?
% ensimmï¿½isenï¿?on yksilï¿? jonka poistaminen korista pop
% parantaisi korin logml:ï¿½ï¿½ eniten, jne...

global COUNTS;      global SUMCOUNTS;

ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind = inds(i);
    rows = globalRows(ind,1):globalRows(ind,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop) = COUNTS(:,:,pop)-diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)-diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop);
    COUNTS(:,:,pop) = COUNTS(:,:,pop)+diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)+diffInSumCounts;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);

%--------------------------------------------------------------------------


function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, globalRows, data)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikï¿?olisi
% muutos logml:ssï¿? mikï¿½li yksilï¿?ind siirretï¿½ï¿½n koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lisï¿½ttï¿½vï¿?
% COUNTS:in siivuun i2, mikï¿½li muutos toteutetaan.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

npops = size(COUNTS,3);
muutokset = LOGDIFF(ind,:);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);
muutokset(i1) = 0;

rows = globalRows(ind,1):globalRows(ind,2);
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = find(muutokset==-Inf);     % Etsitï¿½ï¿½n populaatiot jotka muuttuneet viime kerran jï¿½lkeen.
i2 = setdiff(i2,i1);            
i2_logml = POP_LOGML(i2);

ni2 = length(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[ni2 1]);
new_i2_logml = computePopulationLogml(i2);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[ni2 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;
LOGDIFF(ind,:) = muutokset;


%----------------------------------------------------------------------


function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukumï¿½ï¿½rï¿½t (vastaavasti kuin COUNTS:issa), jotka ovat data:n 
% riveillï¿?rows. rows pitï¿½ï¿½ olla vaakavektori.

diffInCounts = zeros(max_noalle, nloci);
for i=rows
    row = data(i,:);
    notEmpty = find(row>=0);
    
    if length(notEmpty)>0
        diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
            diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
    end
end    

%------------------------------------------------------------------------


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, diffInCounts)
% Suorittaa globaalien muuttujien muutokset, kun yksilï¿?ind
% on siirretï¿½ï¿½n koriin i2.

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
global POP_LOGML;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2]);


%--------------------------------------------------------------------------
%--

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2( ...
    i1, globalRows, data);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikï¿?olisi
% muutos logml:ssï¿? mikï¿½li korin i1 kaikki yksilï¿½t siirretï¿½ï¿½n
% koriin i. 

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
muutokset = zeros(npops,1);

i1_logml = POP_LOGML(i1);

inds = find(PARTITION==i1);
ninds = length(inds);

if ninds==0
    diffInCounts = zeros(size(COUNTS,1), size(COUNTS,2));
    return;
end

rows = [];
for ind = inds
    lisa = globalRows(ind,1):globalRows(ind,2);
    rows = [rows; lisa'];
    %rows = [rows; globalRows{ind}'];
end

diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%---------------------------------------------------------------------------------


function updateGlobalVariables2( ...
    i1, i2, diffInCounts);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilï¿½t siirretï¿½ï¿½n koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;

inds = find(PARTITION==i1);
PARTITION(inds) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2);


%--------------------------------------------------------------------------
%----

function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
    data, i1)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mikï¿?olisi muutos logml:ssï¿? jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirretï¿½ï¿½n koriin j.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
npops2 = length(unique(T2));
muutokset = zeros(npops2, npops);

i1_logml = POP_LOGML(i1);

for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
        rows = [];
        for ind = inds
            lisa = globalRows(ind,1):globalRows(ind,2);
            rows = [rows; lisa'];
            %rows = [rows; globalRows{ind}'];
        end
        diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
        diffInSumCounts = sum(diffInCounts);

        COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
        new_i1_logml = computePopulationLogml(i1);
        COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';
    
        COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
        new_i2_logml = computePopulationLogml(i2)';
        COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

        muutokset(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end    
end



%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, i2);
% Suorittaa globaalien muuttujien pï¿½ivitykset, kun yksilï¿½t 'muuttuvat'
% siirretï¿½ï¿½n koriin i2. Ennen siirtoa yksilï¿½iden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2]);


%----------------------------------------------------------------------------


function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sisï¿½ltï¿½ï¿½ yksilï¿½iden inds2
% vï¿½liset etï¿½isyydet. ninds=kaikkien yksilï¿½iden lukumï¿½ï¿½rï¿?

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

%--------------------------------------------------------------------------


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


%-----------------------------------------------------------------------------------


function npops = poistaTyhjatPopulaatiot(prioriPopLkm)
% Poistaa tyhjentyneet populaatiot COUNTS:ista ja 
% SUMCOUNTS:ista, ADJPRIOR:ista ja SUMPRIOR:ista. 
% Pï¿½ivittï¿½ï¿½ npops:in ja PARTITION:in.

global COUNTS;
global SUMCOUNTS;
global PARTITION;
global ADJPRIOR;
global SUMPRIOR;
global LOGDIFF;

notEmpty = union(find(any(SUMCOUNTS,2)) , 1:prioriPopLkm);
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);

ADJPRIOR = ADJPRIOR(:,:,notEmpty);
SUMPRIOR = SUMPRIOR(notEmpty,:);
LOGDIFF = LOGDIFF(:,notEmpty);

for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
end
npops = length(notEmpty);


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


%-------------------------------------------------------------------------


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

%------------------------------------------------------------------

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

%---------------------------------------------------------------


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


%-------------------------------------------------------------------


function writeTrainedMixtureInfo(logml, rows, data, outPutFile, ...
    inputFile, suIndNames, suNames, pIndNames, pNames, partitionSummary)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global ADJPRIOR;
ninds = size(rows,1);
npops =  size(COUNTS,3);
n_clust_with_su = length(unique(PARTITION));

if length(outPutFile)>0
    fid = fopen(outPutFile,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
end

dispLine;
disp('RESULTS OF TRAINED MIXTURE ANALYSIS:');
disp(['Data file: ' inputFile]);
disp(['Number of clustered groups: ' ownNum2Str(ninds)]);
disp(['Number of populations having prior information: ' ownNum2Str(size(pNames,1))]);
disp(['In the optimal partition the samling units were in ' ownNum2Str(n_clust_with_su) ' clusters.']);
disp(['Log(marginal likelihood) of the optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
    fprintf(fid,'%s \n', ['    ']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['RESULTS OF TRAINED MIXTURE ANALYSIS:']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clustered groups: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of populations having prior information: ' ownNum2Str(size(pNames,1))]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['In the optimal partition the sampling units were in ' ownNum2Str(n_clust_with_su) ' clusters.']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Log(marginal likelihood) of the optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
    fprintf(fid,'\n');
end

%cluster_count = length(unique(PARTITION));
cluster_count = size(COUNTS,3);
disp(['Best Partition: ']);
if (fid ~= -1)
    fprintf(fid,'%s \n',['Best Partition: ']); fprintf(fid,'\n');
end
for m=1:cluster_count
    susInM = find(PARTITION==m);  %Sampling units in pop m.
    text = ['Cluster ' num2str(m) ': {'];
    length_of_beginning = 11 + floor(log10(m));
    
    if m < size(pNames,1)
        % populaatiolle on allokoitu prioriyksilï¿½itï¿?
        text = [text '['];
        k = pNames{m,2};
        text = [text pIndNames{k}{1}];
        for k = pNames{m,2}+1:pNames{m+1,2}-1
            text = [text ', ' pIndNames{k}{1}];
        end
        text = [text '], '];
    elseif m == size(pNames,1)
        text = [text '['];
        k = pNames{m,2};
        text = [text pIndNames{k}{1}];
        for k = pNames{m,2}+1:length(pIndNames)
            text = [text ', ' pIndNames{k}{1}];
        end
        text = [text '], '];
    end       
    
    cluster_size = length(susInM);
    
    for k = 1:cluster_size  % Kï¿½y lï¿½pi m:ï¿½ï¿½n kuuluvat samling unit:it
        text = [text '['];
        su = susInM(k);  % sampling unit su kuuluu populaatioon m.
        ekaNimi = suNames{su,2};
        if su<size(suNames,1)
            vikaNimi = suNames{su+1,2}-1;
        else
            vikaNimi = length(suIndNames);
        end
        
        for ind = ekaNimi : vikaNimi
            if ind==ekaNimi text = [text suIndNames{ind}{1}];
            else text = [text ', ' suIndNames{ind}{1}];
            end
        end
        text = [text '], '];
    end
    text = text(1:end-2);  %Ota pois viimeinen pilkku.
     
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
    disp('Changes in log(marginal likelihood) if sampling unit i is moved to cluster j:');
    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if sampling unit i is moved to cluster j:']); fprintf(fid, '\n');
    end

    ekarivi = 'group      ';
    for i = 1:cluster_count
        ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
    end
    disp(ekarivi);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [ekarivi]); fprintf(fid, '\n');
    end

    for ind = 1:ninds
        [muutokset, diffInCounts] = laskeMuutokset(ind, rows, data);
        rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
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
    for pop1 = 1:npops
        prior = ADJPRIOR(:,:,pop1);
        prior(find(prior==1))=0;
        nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yhtï¿?alleelia.
        prior(1,nollia)=1;
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
%             fprintf(fid, '%s \n', [rivi]); %fprintf(fid, '\n');
%         end
    end

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
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values:');

if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values:']); fprintf(fid, '\n');
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

disp(' ');
disp(' ');
disp('Probabilities for number of clusters: (#clusters: prob)');

if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['Probabilities for number of clusters: (#clusters: prob)']); fprintf(fid, '\n');
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
        line = [num2str(npopsTaulu(i)) ':   ' num2str(probs(i))];
        disp(line);
        if (fid ~= -1)
            fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
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

%--------------------------------------------------------------

function num2 = omaRound(num)
% Pyï¿½ristï¿½ï¿½ luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------

function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:inï¿?
% yks tï¿½ytyy olla kokonaisluku, joka on 
% vï¿½hintï¿½ï¿½n -1:n suuruinen. Pienemmillï¿?
% luvuilla tapahtuu jokin pyï¿½ristysvirhe.

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


%-----------------------------------------------


function ninds = testaaOnkoKunnollinenBapsData(data)
%Tarkastaa onko viimeisessï¿?sarakkeessa kaikki
%luvut 1,2,...,n johonkin n:ï¿½ï¿½n asti.
%Tarkastaa lisï¿½ksi, ettï¿?on vï¿½hintï¿½ï¿½n 2 saraketta.
if size(data,1)<2
    ninds = 0; return;
end
lastCol = data(:,end);
ninds = max(lastCol);
if ~isequal((1:ninds)',unique(lastCol))
    ninds = 0; return;
end



%--------------------------------------------------------------------------
function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(raw_data)
% Alkuperäisen datan viimeinen sarake kertoo, milt?yksilölt?
% kyseinen rivi on peräisin. Funktio tutkii ensin, ett?montako
% rivi?maksimissaan on peräisin yhdelt?yksilölt? jolloin saadaan
% tietää onko kyseess?haploidi, diploidi jne... Tämän jälkeen funktio
% lisää tyhji?rivej?niille yksilöille, joilta on peräisin vähemmän
% rivej?kuin maksimimäär?
%   Mikäli jonkin alleelin koodi on =0, funktio muuttaa tämän alleelin
% koodi pienimmäksi koodiksi, joka isompi kuin mikään käytöss?oleva koodi.
% Tämän jälkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
% koodit saavat arvoja välill?1,...,noalle(j).

% added by Lu Cheng, without modification, 16.02.2010

data = raw_data;
nloci=size(raw_data,2)-1;

dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
    isoinAlleeli = max(max(dataApu));
    dataApu(nollat) = isoinAlleeli+1;
    data(:,1:nloci) = dataApu;
end
% dataApu = []; 
% nollat = []; 
% isoinAlleeli = [];

noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
   %alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(logical(alleelitLokuksessaI>=0));
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
        % data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
        data(logical(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;

nind = max(data(:,end));
nrows = size(data,1);
ncols = size(data,2);
rowsFromInd = zeros(nind,1);
for i=1:nind
    rowsFromInd(i) = length(find(data(:,end)==i));
end
maxRowsFromInd = max(rowsFromInd);
a = -999;
emptyRow = repmat(a, 1, ncols);
lessThanMax = find(rowsFromInd < maxRowsFromInd);
missingRows = maxRowsFromInd*nind - nrows;
data = [data; zeros(missingRows, ncols)];
pointer = 1;
for ind=lessThanMax'    %Käy läpi ne yksilöt, joilta puuttuu rivej?
    miss = maxRowsFromInd-rowsFromInd(ind);  % Tält?yksilölt?puuttuvien lkm.
    for j=1:miss
        rowToBeAdded = emptyRow;
        rowToBeAdded(end) = ind;
        data(nrows+pointer, :) = rowToBeAdded;
        pointer = pointer+1;
    end
end
data = sortrows(data, ncols);   % Sorttaa yksilöiden mukaisesti
newData = data;
rowsFromInd = maxRowsFromInd;

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
end


%-------------------------------------------------------------------------

function [Z, dist] = newGetDistances(data, rowsFromInd)

ninds = max(data(:,end));
nloci = size(data,2)-1;
riviLkm = nchoosek(double(ninds),2);

% empties = find(data<0);
% data(empties)=0;
data(logical(data<0)) = 0;
data = uint16(data);

pariTaulu = zeros(riviLkm,2);
aPointer=1;

for a=1:ninds-1
    pariTaulu(aPointer:aPointer+double(ninds-1-a),1) = ones(ninds-a,1,'uint16')*a;
    pariTaulu(aPointer:aPointer+double(ninds-1-a),2) = uint16((a+1:ninds)');
    aPointer = aPointer+double(ninds-a);
end

eka = pariTaulu(:,ones(1,rowsFromInd));
eka = eka * rowsFromInd;
miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
eka = eka - miinus;

toka = pariTaulu(:,ones(1,rowsFromInd)*2);
toka = toka * rowsFromInd;
toka = toka - miinus;

eka = uint16(eka);
toka = uint16(toka);

clear pariTaulu; clear miinus;

summa = uint16(zeros(riviLkm,1));
vertailuja = uint16(zeros(riviLkm,1));

x = zeros(size(eka));    x = uint16(x);
y = zeros(size(toka));   y = uint16(y);
% fprintf(1,'%%10');
for j=1:nloci;
    
    for k=1:rowsFromInd
        x(:,k) = data(eka(:,k),j);
        y(:,k) = data(toka(:,k),j);
    end

    for a=1:rowsFromInd
        for b=1:rowsFromInd
            vertailutNyt = uint16(x(:,a)>0 & y(:,b)>0);
            vertailuja = vertailuja + vertailutNyt;
            lisays = (x(:,a)~=y(:,b) & vertailutNyt);
            summa = summa + uint16(lisays);
        end
    end
    % fprintf(1,'\b\b');
    % fprintf(1,'%d',floor(10+80*j/nloci));
end

clear x;    clear y;   clear vertailutNyt;
clear eka; clear toka; clear data; clear lisays;
dist = zeros(length(vertailuja),1);
% nollat = find(vertailuja==0);
% dist(nollat) = 1;
dist(logical(vertailuja==0)) = 1;
muut = find(vertailuja>0);
dist(muut) = double(summa(muut))./double(vertailuja(muut));
clear summa; clear vertailuja; clear muut;

Z = computeLinkage(dist');