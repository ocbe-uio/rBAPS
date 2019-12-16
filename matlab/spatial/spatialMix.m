function [logml, npops, partitionSummary]=spatialMix(c,npopsTaulu)
% Greedy search algorithm with unknown number of classes for spatial
% clustering.

logml = 1;

global PARTITION; global COUNTS;
global SUMCOUNTS; 
global SEPCOUNTS; global CLIQCOUNTS;
global LOGDIFF;
clearGlobalVars;

data = c.data; rowsFromInd = c.rowsFromInd; 
noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
cliques = c.cliques; separators = c.separators; rows = c.rows;

if isfield(c, 'dist')
    dist = c.dist; Z = c.Z;
end

clear c;


if nargin < 2;
    npopstext = [];
    ready = false;
    teksti = 'Input upper bound to the number of populations (possibly multiple values): ';
    while ready == false
        npopstextExtra = inputdlg(teksti ,...
            'Input maximum number of populations',1,{'20'});
        drawnow
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
        ykkoset = find(npopsTaulu==1);
        npopsTaulu(ykkoset) = [];   % Mik‰li ykkˆsi‰ annettu yl‰rajaksi, ne poistetaan.
        if isempty(npopsTaulu)
            return
        end
        clear ykkoset;
    end
end

nruns = length(npopsTaulu);

maxnpops = max(npopsTaulu);
logmlBest = -1e50;
partitionSummary = -1e50*ones(30,2,maxnpops);  % Tiedot 30 parhaasta partitiosta (npops ja logml)
partitionSummary(:,1,:) = zeros(30,1,maxnpops);
worstLogml = -1e50*ones(1, maxnpops); worstIndex = ones(1, maxnpops);


initData = data;
data = initData(:, 1:end-1);

for run = 1:nruns 
    npops = npopsTaulu(run);
    dispLine;
    disp(['Run ' num2str(run) '/' num2str(nruns) ...
        ', maximum number of populations ' num2str(npops) '.']);       
    disp(' ');
    disp('Computing initial partition');
    
    [initialPartition, counts, sumcounts] = initSpatialMultiMixture(initData, ...
        npops, Z, rows, noalle, dist, adjprior, priorTerm,0);
    
    PARTITION = initialPartition;    
    [cliqcounts, sepcounts] = computeCounts(cliques, separators, npops);
    
    COUNTS = zeros(max(noalle), size(data,2),npops); SUMCOUNTS = zeros(npops,size(data,2));
    COUNTS(:,:,1:size(counts,3)) = counts; SUMCOUNTS(1:size(sumcounts,1),:) = sumcounts;
    CLIQCOUNTS = cliqcounts; SEPCOUNTS = sepcounts;
    
    ninds = length(PARTITION);
    %maxsize=max([max(noalle) npops]);
    %initializeGammaln(ninds, rowsFromInd, maxsize);
    logml = computeLogml(adjprior, priorTerm);
    
    LOGDIFF = repmat(-Inf,ninds,npops);
    
    nnotEmptyPops = length(unique(PARTITION));
    
    if logml>worstLogml(nnotEmptyPops);
        [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
            partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
        if (added==1)
            [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                min(partitionSummary(:,2,nnotEmptyPops));  
        end
    end
    
    clear initialPartition; clear counts; clear sumcounts; 
    clear cliqcounts; clear sepcounts;

    % PARHAAN MIXTURE-PARTITION ETSIMINEN
    roundTypes = [1];  %Ykkˆsvaiheen sykli kahteen kertaan.
    nRoundTypes = 7;
    kokeiltu = zeros(nRoundTypes, 1);
    ready = 0; vaihe = 1;
    
    
    disp(['Mixture analysis started with initial ' num2str(nnotEmptyPops) ' populations.']);

    while ready ~= 1
	    muutoksia = 0;   
	
        disp(['Performing steps: ' num2str(roundTypes)]);
    
        for n = 1:length(roundTypes)
        
            round = roundTypes(n);
                        
            if kokeiltu(round) == 1     %Askelta kokeiltu viime muutoksen j‰lkeen
                                
            elseif round==0 | round==1   %Yksilˆn siirt‰minen toiseen populaatioon.
                
                inds = randperm(ninds);
                muutosNyt = 0;
                for ind = inds
                    i1 = PARTITION(ind);
                    [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
                        data, adjprior, priorTerm, logml, cliques, separators);
            
                    if round==1, [maxMuutos, i2] = max(muutokset);
                    elseif round==0, [maxMuutos, i2] = arvoSeuraavaTila(muutokset, logml);
                    end
            
                    if (i1~=i2 & maxMuutos>1e-5)
                        % Tapahtui muutos
                        muutoksia = 1;
                        if muutosNyt == 0
                            disp('action 1');
                            muutosNyt = 1;
                            kokeiltu = zeros(nRoundTypes,1);
                        end
                        updateGlobalVariables(ind, i2, diffInCounts,...
                            cliques, separators, adjprior, priorTerm);
                        logml = logml+maxMuutos;
                                            
                        nnotEmptyPops = length(unique(PARTITION));

                        if logml>worstLogml(nnotEmptyPops);
                            [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
                                partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
                            if (added==1)
                                [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                                    min(partitionSummary(:,2,nnotEmptyPops)); 
                            end
                        end
                    end
                end
                if muutosNyt == 0
                   kokeiltu(round) = 1;
                end 
                
                                                
            elseif round==2  %Populaation yhdist‰minen toiseen.
                maxMuutos = 0;
                
                for pop = 1:npops
                    [muutokset, diffInCounts] = laskeMuutokset2(pop, rows, ...
                        data, adjprior, priorTerm, logml, cliques, separators);
                    [isoin, indeksi] = max(muutokset);
                    if isoin>maxMuutos
                        maxMuutos = isoin;
                        i1 = pop;
                        i2 = indeksi;
                        diffInCountsBest = diffInCounts;
                    end
                end
                
                
                
                if maxMuutos > 1e-5
                    disp('action 2');
                    muutoksia = 1;
                    kokeiltu = zeros(nRoundTypes,1);
                    updateGlobalVariables2(i1,i2, diffInCountsBest, ...
                        cliques, separators, adjprior, priorTerm);
                    logml=logml + maxMuutos;
                    
                    nnotEmptyPops = length(unique(PARTITION));
    
                    if logml>worstLogml(nnotEmptyPops);
                        [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
                            partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
                        if (added==1)
                            [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                                min(partitionSummary(:,2,nnotEmptyPops)); 
                        end
                    end
                else
                    kokeiltu(round) = 1;
                end
                
            
            elseif round==3 | round==4 %Populaation jakaminen osiin.
                maxMuutos = 0;
                for pop = 1:npops
                    inds2 = find(PARTITION==pop);
                    ninds2 = length(inds2);
                    if ninds2>5
                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = linkage(dist2');
                        if round==3
                            npops2 = min(20, floor(ninds2 / 5));  %Moneenko osaan jaetaan
                        elseif round==4
                            npops2 = 2;
                        end
                        T2 = cluster_own(Z2, npops2);
                        
                        muutokset = laskeMuutokset3(T2, inds2, rows, data, ...
                            adjprior, priorTerm, pop, logml, cliques, separators);
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
                    if round==3
                        disp('action 3');
                    else
                        disp('action 4');
                    end
                    kokeiltu = zeros(nRoundTypes,1);
                    rivit = [];
                    for i = 1:length(muuttuvat)
                        ind = muuttuvat(i);
                        lisa = rows(ind,1):rows(ind,2);
                        rivit = [rivit lisa];
                    end
                    diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
                        size(COUNTS,2), data);
                    i1 = PARTITION(muuttuvat(1));
                    updateGlobalVariables3(muuttuvat, diffInCounts, ...
                        adjprior, priorTerm, i2, cliques, separators);
                    logml = logml + maxMuutos;
                    
                    nnotEmptyPops = length(unique(PARTITION));

                    if logml>worstLogml(nnotEmptyPops);
                        [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
                            partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
                        if (added==1)
                            [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                                min(partitionSummary(:,2,nnotEmptyPops));  end
                    end
                else
                    kokeiltu(round)=1;
                end

            elseif round == 5 | round == 6
                j=0;
                muutettu = 0;
                %poplogml = POP_LOGML;
                partition = PARTITION;
                counts = COUNTS;
                sumcounts = SUMCOUNTS;
                cliqcounts = CLIQCOUNTS;
                sepcounts = SEPCOUNTS;
                logdiff = LOGDIFF;
                
                pops = randperm(npops);            
                while (j < npops & muutettu == 0)
                    j = j+1;
                    pop = pops(j);
                    totalMuutos = 0;
                    inds = find(PARTITION==pop);
                    if round == 5
                        aputaulu = [inds rand(length(inds),1)];
                        aputaulu = sortrows(aputaulu,2);
                        inds = aputaulu(:,1)';
                    elseif round == 6
                        inds = returnInOrder(inds, pop, ...
                            rows, data, adjprior, priorTerm);
                    end
                    
                    i=0;
                                
                    while (length(inds) > 0 & i < length(inds))
                        i = i+1;
                        ind =inds(i);                    
                                                 
                        [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
                            data, adjprior, priorTerm, logml, cliques, separators);
                        muutokset(pop) = -1e50;   % Varmasti ei suurin!!!
                        [maxMuutos, i2] = max(muutokset);
                                                
                        updateGlobalVariables(ind, i2, diffInCounts, ...
                            cliques, separators, adjprior, priorTerm);
                        
                        totalMuutos = totalMuutos+maxMuutos;
                        logml = logml+maxMuutos;
                        if round == 6
                            % Lopetetaan heti kun muutos on positiivinen.
                            if totalMuutos > 1e-5
                                i=length(inds);
                            end
                        end
                    end
                
                   if totalMuutos>1e-5
                        if round == 5
                            disp('action 5');
                        elseif round == 6
                            disp('action 6');
                        end
                        kokeiltu = zeros(nRoundTypes,1);
                        muutettu=1;
                        muutoksia = 1;  % Ulompi kirjanpito.
                        nnotEmptyPops = length(unique(PARTITION));

                        if logml>worstLogml(nnotEmptyPops);
                            [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
                                partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
                            if (added==1)
                                [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                                    min(partitionSummary(:,2,nnotEmptyPops));
                            end
                        end
                   else
                        % Miss‰‰n vaiheessa tila ei parantunut.
                        % Perutaan kaikki muutokset.
                        PARTITION = partition;
                        SUMCOUNTS = sumcounts;
                        %POP_LOGML = poplogml;
                        COUNTS = counts;
                        logml = logml - totalMuutos;
                        CLIQCOUNTS = cliqcounts;
                        SEPCOUNTS = sepcounts;
                        LOGDIFF = logdiff;
                        kokeiltu(round)=1;
                    end
                end
                clear partition; clear sumcounts; clear counts; 
                clear cliqcounts; clear sepcounts;
                
            elseif round == 7
                emptyPop = findEmptyPop(npops);
                j = 0;
                pops = randperm(npops);
                muutoksiaNyt = 0;
                if emptyPop == -1
                    j = npops;
                end
                while (j < npops)
                    j = j +1;
                    pop = pops(j);
                    inds2 = find(PARTITION==pop);
                    ninds2 = length(inds2);
                    if ninds2 > 5
                        partition = PARTITION;
                        sumcounts = SUMCOUNTS;
                        counts = COUNTS;
                        cliqcounts = CLIQCOUNTS;
                        sepcounts = SEPCOUNTS;
                        oldLogml = logml;
                        logdiff = LOGDIFF;
                        
                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = linkage(dist2');
                        T2 = cluster_own(Z2, 2);
                        muuttuvat = inds2(find(T2==1));

                        rivit = [];
                        for i = 1:length(muuttuvat)
                            ind = muuttuvat(i);
                            lisa = rows(ind,1):rows(ind,2);
                            rivit = [rivit lisa];
                        end
                        diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
                            size(COUNTS,2), data);
                        updateGlobalVariables3(muuttuvat, diffInCounts, ...
                            adjprior, priorTerm, emptyPop, cliques, separators);
                        
                        logml = computeLogml(adjprior, priorTerm);
                        
                        muutettu = 1;
                        while (muutettu == 1)
                            muutettu = 0;
                            % Siirret‰‰n yksilˆit‰ populaatioiden v‰lill‰
                            muutokset = laskeMuutokset5(inds2, rows, data, ...
                                adjprior, priorTerm, logml, cliques, separators, pop, emptyPop);
                            
                            [maxMuutos, indeksi] = max(muutokset);                                                       
                            muuttuva = inds2(indeksi);
                            
                            if (PARTITION(muuttuva) == pop)
                                i2 = emptyPop;
                            else
                                i2 = pop;
                            end
                                                     
                            if maxMuutos > 1e-5
                                rivit = rows(muuttuva,1):rows(muuttuva,2);
                                diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
                                    size(COUNTS,2), data);
                                updateGlobalVariables3(muuttuva, diffInCounts, ...
                                    adjprior, priorTerm, i2, cliques, separators);
                                muutettu = 1;
                                logml = logml + maxMuutos;
                            end
                            
                        end
                        
                        if logml > oldLogml + 1e-5
                            muutoksia = 1;
                            nnotEmptyPops = length(unique(PARTITION));

                            if logml>worstLogml(nnotEmptyPops);
                                [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
                                    partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
                                if (added==1)
                                    [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
                                        min(partitionSummary(:,2,nnotEmptyPops));  end
                            end
                            if muutoksiaNyt == 0
                                disp('action 7');
                                muutoksiaNyt = 1;
                            end
                            kokeiltu = zeros(nRoundTypes,1);
                            j = npops;  % Lopetetaan, koska muutos syntyi
                        else
                            %palutetaan vanhat arvot
                            PARTITION = partition;
                            SUMCOUNTS = sumcounts;
                            COUNTS = counts;
                            CLIQCOUNTS = cliqcounts;
                            SEPCOUNTS = sepcounts;
                            LOGDIFF = logdiff;
                            logml = oldLogml;
                        end
                            
                    end
                    
                end
                
                if muutoksiaNyt == 0
                    kokeiltu(round)=1;
                end
                    
            end
        end
    
    
        if muutoksia == 0
            if vaihe==1
                vaihe = 2;
            elseif vaihe==2
                vaihe = 3;
            elseif vaihe==3
                vaihe = 4;
            elseif vaihe==4;
                vaihe = 5;
            elseif vaihe==5
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
                roundTypes=[5 5 7];
            elseif vaihe==4
                roundTypes=[4 3 1];
            elseif vaihe==5
                roundTypes=[6 2 7 3 4 1];
            end
        end
    end

    % TALLENNETAAN
    
    npops = poistaTyhjatPopulaatiot(npops);
    %POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    disp(['Found partition with ' num2str(npops) ' populations.']);
    disp(['Log(ml) = ' num2str(logml)]);
    disp(' ');
    
    if logml>logmlBest
        % P‰ivitet‰‰n parasta lˆydetty‰ partitiota.
        logmlBest = logml;
        npopsBest = npops;
        partitionBest = PARTITION;
        countsBest = COUNTS;
        sumCountsBest = SUMCOUNTS;
        %pop_logmlBest = POP_LOGML;
        cliqCountsBest = CLIQCOUNTS;
        sepCountsBest = SEPCOUNTS;
        logdiffbest = LOGDIFF;
    end
end
    
logml = logmlBest;
npops = npopsBest;
PARTITION = partitionBest;
COUNTS = countsBest;
SUMCOUNTS = sumCountsBest;
%POP_LOGML = pop_logmlBest;
CLIQCOUNTS = cliqCountsBest;
SEPCOUNTS = sepCountsBest;
LOGDIFF = logdiffbest;

%-------------------------------------------------------------------------------------

function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
%global POP_LOGML; POP_LOGML = [];
global SEPCOUNTS; SEPCOUNTS = [];
global CLIQCOUNTS; CLIQCOUNTS = [];
global LOGDIFF; LOGDIFF = [];

%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedet‰‰n, ett‰ annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss‰ ei viel‰ ole
% annettua logml arvoa, niin lis‰t‰‰n worstIndex:in kohtaan uusi logml ja
% nykyist‰ partitiota vastaava nclusters:in arvo. Muutoin ei tehd‰ mit‰‰n.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lˆydetty partitio ei ole viel‰ kirjattuna summaryyn.
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
% Suorittaa globaalien muuttujien muutokset, kun yksilˆ ind
% siirret‰‰n koriin i2.

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
global CLIQCOUNTS;
global SEPCOUNTS;
global LOGDIFF;

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

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;

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
global LOGDIFF;

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

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, ...
    adjprior, priorTerm, i2, cliques, separators);
% Suorittaa globaalien muuttujien p‰ivitykset, kun yksilˆt 'muuttuvat'
% siirret‰‰n koriin i2. Ennen siirtoa yksilˆiden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;      global CLIQCOUNTS;
global SUMCOUNTS;   global SEPCOUNTS;
global LOGDIFF;
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

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;

%POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%----------------------------------------------------------------------


function inds = returnInOrder(inds, pop, globalRows, data, ...
    adjprior, priorTerm)
% Palauttaa yksilˆt j‰rjestyksess‰ siten, ett‰ ensimm‰isen‰ on
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
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik‰ olisi
% muutos logml:ss‰, mik‰li yksilˆt inds siirret‰‰n koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lis‰tt‰v‰
% COUNTS:in siivuun i2, mik‰li muutos toteutetaan.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   %global POP_LOGML;
global CLIQCOUNTS;  global SEPCOUNTS;
global LOGDIFF;

npops = size(COUNTS,3);
muutokset = LOGDIFF(ind,:);

counts = COUNTS;
sumcounts = SUMCOUNTS;

[emptyPop, pops] = findEmptyPop(npops);

i1 = PARTITION(ind);
muutokset(i1) = 0;

i2 = [pops(find(pops~=i1))];
if emptyPop > 0
    i2 =[i2 emptyPop];
end

i2 = sort(i2);
laskematta = find(muutokset==-Inf);
i2 = intersect(i2,laskematta);

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
LOGDIFF(ind,:) = muutokset;

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2(i1, globalRows, ...
    data, adjprior, priorTerm, logml, cliques, separators);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik‰ olisi
% muutos logml:ss‰, mik‰li korin i1 kaikki yksilˆt siirret‰‰n
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
% kertoo, mik‰ olisi muutos logml:ss‰, jos populaation i1 osapopulaatio
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

% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik‰ olisi
% muutos logml:ss‰, mik‰li yksilˆ i vaihtaisi koria i1:n ja i2:n v‰lill‰.
    
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
        pop1 = i1;  %mist‰
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
% riveill‰ rows.

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
global LOGDIFF;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);
CLIQCOUNTS = CLIQCOUNTS(:,notEmpty);
SEPCOUNTS = SEPCOUNTS(:,notEmpty);
LOGDIFF = LOGDIFF(:,notEmpty);

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
% Laskee arvot cliqcounts:lle ja sepcounts:lle

function [cliqcounts, sepcounts] = computeCounts(cliques, separators, npops)

global PARTITION;
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
% Laskee muutoksen CLIQCOUNTS:ssa (tai SEPCOUNTS:ssa, jos syˆtteen‰
% separators) kun yksilˆt inds siirret‰‰n.
% diffInCliqcounts on ncliq*1 taulu, joka on CLIQCOUNTS:n sarakkeesta josta
% yksilˆt inds siirret‰‰n ja lis‰tt‰v‰ sarakkeeseen, johon yksilˆt
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

function [logml, spatialPrior] = computeLogml(adjprior,priorTerm)

%global GAMMA_LN;
global CLIQCOUNTS;
global SEPCOUNTS;
global PARTITION;

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


global COUNTS;
global SUMCOUNTS;
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
% v‰liset et‰isyydet. ninds=kaikkien yksilˆiden lukum‰‰r‰.

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
function mjono = logml2String(logml)
% Palauttaa logml:n string-esityksen.

mjono = '       ';
if abs(logml)<10000
    %Ei tarvita e-muotoa
    mjono(7) = palautaYks(abs(logml),-1);
    mjono(6) = '.';
    mjono(5) = palautaYks(abs(logml),0);
    mjono(4) = palautaYks(abs(logml),1);
    mjono(3) = palautaYks(abs(logml),2);
    mjono(2) = palautaYks(abs(logml),3);
    pointer = 2;
    while mjono(pointer)=='0' & pointer<7
        mjono(pointer) = ' ';
        pointer=pointer+1;
    end
    if logml<0
        mjono(pointer-1) = '-';
    end 
else
    suurinYks = 4;
    while abs(logml)/(10^(suurinYks+1)) >= 1
        suurinYks = suurinYks+1;
    end
    if suurinYks<10
        mjono(7) = num2str(suurinYks);
        mjono(6) = 'e';
        mjono(5) = palautaYks(abs(logml),suurinYks-1);
        mjono(4) = '.';
        mjono(3) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(2) = '-';
        end
    elseif suurinYks>=10
        mjono(6:7) = num2str(suurinYks);
        mjono(5) = 'e';
        mjono(4) = palautaYks(abs(logml),suurinYks-1);
        mjono(3) = '.';
        mjono(2) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(1) = '-';
        end
    end
end

function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in‰ 
% yks t‰ytyy olla kokonaisluku, joka on 
% v‰hint‰‰n -1:n suuruinen. Pienemmill‰
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
% Palauttaa ensimm‰isen tyhj‰n populaation indeksin. Jos tyhji‰
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

function [sumcounts, counts, logml] = ...
    initialCounts(partition, data, npops, rowsFromInd, noalle)

nloci=size(data,2);
ninds = size(data,1)/rowsFromInd;

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

%---------------------------------------------------------------


function dispLine;
disp('---------------------------------------------------');
