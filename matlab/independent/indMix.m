function [logml, npops, partitionSummary]=indMix(c,npops,dispText)
% Greedy search algorithm with unknown number of classes for regular
% clustering.
% Input npops is not used if called by greedyMix or greedyPopMix.

logml = 1;

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;
clearGlobalVars;

noalle = c.noalle; rows = c.rows; data = c.data;
adjprior = c.adjprior; priorTerm = c.priorTerm; rowsFromInd = c.rowsFromInd;

if isfield(c,'dist')
    dist = c.dist; Z = c.Z;
end

clear c;

if nargin < 2
    dispText = 1;
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
        npopsTaulu(ykkoset) = [];   % Mikäli ykkösiä annettu ylärajaksi, ne poistetaan.
        if isempty(npopsTaulu)
            logml = 1; partitionSummary=1; npops=1;
            return
        end
        clear ykkoset;
    end
else
    npopsTaulu = npops;
end

nruns = length(npopsTaulu);

initData = data;
data = data(:,1:end-1);

logmlBest = -1e50;
partitionSummary = -1e50*ones(30,2);  % Tiedot 30 parhaasta partitiosta (npops ja logml)
partitionSummary(:,1) = zeros(30,1);
worstLogml = -1e50; worstIndex = 1;

for run = 1:nruns
    npops = npopsTaulu(run);
    if dispText
        dispLine;
        disp(['Run ' num2str(run) '/' num2str(nruns) ...
            ', maximum number of populations ' num2str(npops) '.']);
    end
    ninds = size(rows,1);

    initialPartition = admixture_initialization(initData, npops, Z);
    [sumcounts, counts, logml] = ...
        initialCounts(initialPartition, data, npops, rows, noalle, adjprior);
    PARTITION = zeros(ninds, 1);
    for i=1:ninds
        apu = rows(i);
        PARTITION(i) = initialPartition(apu(1));
    end

    COUNTS = counts; SUMCOUNTS = sumcounts;
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    LOGDIFF = repmat(-Inf,ninds,npops);
    clear initialPartition; clear counts; clear sumcounts;

    % PARHAAN MIXTURE-PARTITION ETSIMINEN
    nRoundTypes = 7;
    kokeiltu = zeros(nRoundTypes, 1);
    roundTypes = [1 1];  %Ykkösvaiheen sykli kahteen kertaan.
    ready = 0; vaihe = 1;

    if dispText
        disp(' ');
        disp(['Mixture analysis started with initial ' num2str(npops) ' populations.']);
    end

    while ready ~= 1
        muutoksia = 0;

        if dispText
            disp(['Performing steps: ' num2str(roundTypes)]);
        end

        for n = 1:length(roundTypes)

            round = roundTypes(n);
            kivaluku=0;

            if kokeiltu(round) == 1     %Askelta kokeiltu viime muutoksen jälkeen

            elseif round==0 | round==1   %Yksilön siirtäminen toiseen populaatioon.
                inds = 1:ninds;
                aputaulu = [inds' rand(ninds,1)];
                aputaulu = sortrows(aputaulu,2);
                inds = aputaulu(:,1)';
                muutosNyt = 0;

                for ind = inds
                    i1 = PARTITION(ind);
                    [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
                        data, adjprior, priorTerm);

                    if round==1, [maxMuutos, i2] = max(muutokset);
                    end

                    if (i1~=i2 & maxMuutos>1e-5)
                        % Tapahtui muutos
                        muutoksia = 1;
                        if muutosNyt == 0
                            muutosNyt = 1;
                            if dispText
                                disp('Action 1');
                            end
                        end
                        kokeiltu = zeros(nRoundTypes,1);
                        kivaluku = kivaluku+1;
                        updateGlobalVariables(ind, i2, diffInCounts,...
                            adjprior, priorTerm);
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

            elseif round==2  %Populaation yhdistäminen toiseen.
                maxMuutos = 0;
                for pop = 1:npops
                    [muutokset, diffInCounts] = laskeMuutokset2(pop, rows, ...
                        data, adjprior, priorTerm);
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
                    kokeiltu = zeros(nRoundTypes,1);
                    updateGlobalVariables2(i1,i2, diffInCountsBest, ...
                        adjprior, priorTerm);
                    logml = logml + maxMuutos;
                    if dispText
                        disp('Action 2');
                    end
                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                        if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                    end
                else
                    kokeiltu(round) = 1;
                end


            elseif round==3 || round==4 %Populaation jakaminen osiin.
                maxMuutos = 0;
                ninds = size(rows,1);
                for pop = 1:npops
                    inds2 = find(PARTITION==pop);
                    ninds2 = length(inds2);
                    if ninds2>2
                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = linkage(dist2');
                        if round==3
                            npops2 = max(min(20, floor(ninds2/5)),2);
                        elseif round==4
                            npops2 = 2;  %Moneenko osaan jaetaan
                        end
                        T2 = cluster_own(Z2, npops2);
                        muutokset = laskeMuutokset3(T2, inds2, rows, data, ...
                            adjprior, priorTerm, pop);
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
                    kokeiltu = zeros(nRoundTypes,1);
                    %rows = computeRows(rowsFromInd, muuttuvat, length(muuttuvat));
                    rivit = [];
                    for i = 1:length(muuttuvat)
                        ind = muuttuvat(i);
                        lisa = rows(ind,1):rows(ind,2);
                        rivit = [rivit; lisa'];
                        %rivit = [rivit; rows(ind)'];
                    end
                    diffInCounts = computeDiffInCounts(rivit', size(COUNTS,1), ...
                        size(COUNTS,2), data);
                    i1 = PARTITION(muuttuvat(1));
                    updateGlobalVariables3(muuttuvat, diffInCounts, ...
                        adjprior, priorTerm, i2);
                    logml = logml + maxMuutos;
                    if dispText
                        if round==3
                            disp('Action 3');
                        else
                            disp('Action 4');
                        end
                    end
                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                        if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                    end

                else
                    kokeiltu(round)=1;
                end
            elseif round == 5 || round == 6
                j=0;
                muutettu = 0;
                poplogml = POP_LOGML;
                partition = PARTITION;
                counts = COUNTS;
                sumcounts = SUMCOUNTS;
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

                    i = 0;

                    while (length(inds) > 0 & i < length(inds))
                        i = i+1;
                        ind =inds(i);

                        [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
                            data, adjprior, priorTerm);
                        muutokset(pop) = -1e50;   % Varmasti ei suurin!!!
                        [maxMuutos, i2] = max(muutokset);
                        updateGlobalVariables(ind, i2, diffInCounts,...
                            adjprior, priorTerm);

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
                        kokeiltu = zeros(nRoundTypes,1);
                        muutettu=1;
                        if muutoksia==0
                            muutoksia = 1;  % Ulompi kirjanpito.
                            if dispText
                                if round==5
                                    disp('Action 5');
                                else
                                    disp('Action 6');
                                end
                            end
                        end
                        if logml>worstLogml
                            [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                            if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                        end
                    else
                        % Missään vaiheessa tila ei parantunut.
                        % Perutaan kaikki muutokset.
                        PARTITION = partition;
                        SUMCOUNTS = sumcounts;
                        POP_LOGML = poplogml;
                        COUNTS = counts;
                        logml = logml - totalMuutos;
                        LOGDIFF = logdiff;
                        kokeiltu(round)=1;
                    end
                end
                clear partition; clear sumcounts; clear counts; clear poplogml;

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
                    inds2 = find(PARTITION == pop);
                    ninds2 = length(inds2);
                    if ninds2 > 5
                        partition = PARTITION;
                        sumcounts = SUMCOUNTS;
                        counts = COUNTS;
                        poplogml = POP_LOGML;
                        logdiff = LOGDIFF;

                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = linkage(dist2');
                        T2 = cluster_own(Z2, 2);
                        muuttuvat = inds2(find(T2 == 1));

                        muutokset = laskeMuutokset3(T2, inds2, rows, data, ...
                            adjprior, priorTerm, pop);
                        totalMuutos = muutokset(1, emptyPop);

                        rivit = [];
                        for i = 1:length(muuttuvat)
                            ind = muuttuvat(i);
                            lisa = rows(ind,1):rows(ind,2);
                            rivit = [rivit lisa];
                        end
                        diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
                            size(COUNTS,2), data);

                        updateGlobalVariables3(muuttuvat, diffInCounts, ...
                            adjprior, priorTerm, emptyPop);

                        muutettu = 1;
                        while (muutettu == 1)
                            muutettu = 0;
                            % Siirretään yksilöitä populaatioiden välillä
                            muutokset = laskeMuutokset5(inds2, rows, data, ...
                                adjprior, priorTerm, pop, emptyPop);

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
                                updateGlobalVariables3(muuttuva,diffInCounts, ...
                                    adjprior, priorTerm, i2);
                                muutettu = 1;
                                totalMuutos = totalMuutos + maxMuutos;
                            end

                        end

                        if totalMuutos > 1e-5
                            muutoksia = 1;
                            logml = logml + totalMuutos;
                            if logml>worstLogml
                                [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                                if (added==1)  [worstLogml, worstIndex] = min(partitionSummary(:,2));  end
                            end
                            if muutoksiaNyt == 0
                                if dispText
                                    disp('Action 7');
                                end
                                muutoksiaNyt = 1;
                            end
                            kokeiltu = zeros(nRoundTypes,1);
                            j = npops;
                        else
                            %palutetaan vanhat arvot
                            PARTITION = partition;
                            SUMCOUNTS = sumcounts;
                            COUNTS = counts;
                            POP_LOGML = poplogml;
                            LOGDIFF = logdiff;
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
                roundTypes = [2 1];
            elseif vaihe==3
                roundTypes=[5 5 7];
            elseif vaihe==4
                roundTypes=[4 3 1];
            elseif vaihe==5
                roundTypes=[6 7 2 3 4 1];
            end
        end
    end

    % TALLENNETAAN

    npops = poistaTyhjatPopulaatiot(npops);
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    if dispText
        disp(['Found partition with ' num2str(npops) ' populations.']);
        disp(['Log(ml) = ' num2str(logml)]);
        disp(' ');
    end

    if logml>logmlBest
        % Päivitetään parasta löydettyä partitiota.
        logmlBest = logml;
        npopsBest = npops;
        partitionBest = PARTITION;
        countsBest = COUNTS;
        sumCountsBest = SUMCOUNTS;
        pop_logmlBest = POP_LOGML;
        logdiffbest = LOGDIFF;
    end

end

logml = logmlBest;
npops = npopsBest;
PARTITION = partitionBest;
COUNTS = countsBest;
SUMCOUNTS = sumCountsBest;
POP_LOGML = pop_logmlBest;
LOGDIFF = logdiffbest;

%--------------------------------------------------------------------------

function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];
global LOGDIFF; LOGDIFF = [];

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


function [sumcounts, counts, logml] = ...
    initialPopCounts(data, npops, rows, noalle, adjprior)

nloci=size(data,2);
counts = zeros(max(noalle),nloci,npops);
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

logml = laskeLoggis(counts, sumcounts, adjprior);

%-----------------------------------------------------------------------


function loggis = laskeLoggis(counts, sumcounts, adjprior)
npops = size(counts,3);

logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior))) - ...
    sum(sum(gammaln(1+sumcounts)));
loggis = logml2;


%------------------------------------------------------------------------------------


function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on määritelty pops-muuttujalla.

global COUNTS;
global SUMCOUNTS;
x = size(COUNTS,1);
y = size(COUNTS,2);
z = length(pops);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;
%--------------------------------------------------------------------------


function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, globalRows, data, adjprior, priorTerm)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikä olisi
% muutos logml:ssä, mikäli yksilö ind siirretään koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lisättävä
% COUNTS:in siivuun i2, mikäli muutos toteutetaan.
%
% Lisäys 25.9.2007:
% Otettu käyttöön globaali muuttuja LOGDIFF, johon on tallennettu muutokset
% logml:ssä siirrettäessä yksilöitä toisiin populaatioihin.

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
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = find(muutokset==-Inf);     % Etsitään populaatiot jotka muuttuneet viime kerran jälkeen.
i2 = setdiff(i2,i1);
i2_logml = POP_LOGML(i2);

ni2 = length(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[ni2 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[ni2 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;
LOGDIFF(ind,:) = muutokset;


%----------------------------------------------------------------------


function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukumäärät (vastaavasti kuin COUNTS:issa), jotka ovat data:n
% riveillä rows. rows pitää olla vaakavektori.

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


function updateGlobalVariables(ind, i2, diffInCounts, ...
    adjprior, priorTerm)
% Suorittaa globaalien muuttujien muutokset, kun yksilö ind
% on siirretään koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%--------------------------------------------------------------------------
%--

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2( ...
    i1, globalRows, data, adjprior, priorTerm);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikä olisi
% muutos logml:ssä, mikäli korin i1 kaikki yksilöt siirretään
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
for i = 1:ninds
    ind = inds(i);
    lisa = globalRows(ind,1):globalRows(ind,2);
    rows = [rows; lisa'];
    %rows = [rows; globalRows{ind}'];
end

diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%---------------------------------------------------------------------------------


function updateGlobalVariables2( ...
    i1, i2, diffInCounts, adjprior, priorTerm);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilöt siirretään koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

inds = find(PARTITION==i1);
PARTITION(inds) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2, adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%--------------------------------------------------------------------------
%----

function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
    data, adjprior, priorTerm, i1)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mikä olisi muutos logml:ssä, jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirretään koriin j.

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
        for i = 1:ninds
            ind = inds(i);
            lisa = globalRows(ind,1):globalRows(ind,2);
            rows = [rows; lisa'];
            %rows = [rows; globalRows{ind}'];
        end
        diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
        diffInSumCounts = sum(diffInCounts);

        COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
        new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
        COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';

        COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
        new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
        COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

        muutokset(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end
end

%------------------------------------------------------------------------------------

function muutokset = laskeMuutokset5(inds, globalRows, data, adjprior, ...
    priorTerm, i1, i2)

% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mikä olisi
% muutos logml:ssä, mikäli yksilö i vaihtaisi koria i1:n ja i2:n välillä.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;

ninds = length(inds);
muutokset = zeros(ninds,1);

i1_logml = POP_LOGML(i1);
i2_logml = POP_LOGML(i2);

for i = 1:ninds
    ind = inds(i);
    if PARTITION(ind)==i1
        pop1 = i1;  %mistä
        pop2 = i2;  %mihin
    else
        pop1 = i2;
        pop2 = i1;
    end
    rows = globalRows(ind,1):globalRows(ind,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)-diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)-diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)+diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)+diffInSumCounts;

    new_logmls = computePopulationLogml([i1 i2], adjprior, priorTerm);
    muutokset(i) = sum(new_logmls);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)+diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)+diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)-diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)-diffInSumCounts;
end

muutokset = muutokset - i1_logml - i2_logml;

%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, ...
    adjprior, priorTerm, i2);
% Suorittaa globaalien muuttujien päivitykset, kun yksilöt 'muuttuvat'
% siirretään koriin i2. Ennen siirtoa yksilöiden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%----------------------------------------------------------------------------


function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sisältää yksilöiden inds2
% väliset etäisyydet. ninds=kaikkien yksilöiden lukumäärä.

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


%-----------------------------------------------------------------------------------


function npops = poistaTyhjatPopulaatiot(npops)
% Poistaa tyhjentyneet populaatiot COUNTS:ista ja
% SUMCOUNTS:ista. Päivittää npops:in ja PARTITION:in.

global COUNTS;
global SUMCOUNTS;
global PARTITION;
global LOGDIFF;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);
LOGDIFF = LOGDIFF(:,notEmpty);

for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
end
npops = length(notEmpty);


%---------------------------------------------------------------


function dispLine;
disp('---------------------------------------------------');

%--------------------------------------------------------------

function num2 = omaRound(num)
% Pyöristää luvun num 1 desimaalin tarkkuuteen
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
% string:inä
% yks täytyy olla kokonaisluku, joka on
% vähintään -1:n suuruinen. Pienemmillä
% luvuilla tapahtuu jokin pyöristysvirhe.

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
%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

function initial_partition=admixture_initialization(data_matrix,nclusters, Z)

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

%--------------------------------------------------------------------------

function [sumcounts, counts, logml] = ...
    initialCounts(partition, data, npops, rows, noalle, adjprior)

nloci=size(data,2);
ninds = size(rows, 1);

koot = rows(:,1) - rows(:,2) + 1;
maxSize = max(koot);

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

%initializeGammaln(ninds, maxSize, max(noalle));

logml = laskeLoggis(counts, sumcounts, adjprior);

%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedetään, että annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ssä ei vielä ole
% annettua logml arvoa, niin lisätään worstIndex:in kohtaan uusi logml ja
% nykyistä partitiota vastaava nclusters:in arvo. Muutoin ei tehdä mitään.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt löydetty partitio ei ole vielä kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end

%--------------------------------------------------------------------------

function inds = returnInOrder(inds, pop, globalRows, data, ...
    adjprior, priorTerm)
% Palauttaa yksilöt järjestyksessä siten, että ensimmäisenä on
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

%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimmäisen tyhjän populaation indeksin. Jos tyhjiä
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;
pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end
