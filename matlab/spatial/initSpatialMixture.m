function [partition, counts, sumcounts] = initSpatialMixture(initData, ... 
    npops, Z, rowsFromInd, noalle, dist, adjprior, priorTerm);
% Etsii spatial mixturelle alkutilan baps 3.1:n ahneella algoritmilla.

global PARTITION_IN; global COUNTS_IN;
global SUMCOUNTS_IN; global POP_LOGML_IN;

data = initData(:,1:end-1);
initialPartition = admixture_initialization(initData, npops, Z);

[sumcounts, counts, logml] = ...
    initialCounts(initialPartition, data, npops, rowsFromInd, noalle);

PARTITION_IN = initialPartition(1:rowsFromInd:end);
COUNTS_IN = counts; SUMCOUNTS_IN = sumcounts;

partition = PARTITION_IN;
return

POP_LOGML_IN = computePopulationLogml(1:npops, adjprior, priorTerm);

  
clear initialPartition; clear counts; clear sumcounts; 

% PARHAAN MIXTURE-PARTITION_IN ETSIMINEN
roundTypes = [1 1];  %Ykkösvaiheen sykli kahteen kertaan.
ready = 0; vaihe = 1;
ninds = size(data,1)/rowsFromInd;



while ready ~= 1
    muutoksia = 0;   
	
    for n = 1:length(roundTypes)
        
        round = roundTypes(n);
        kivaluku=0;
        
        if round==0 | round==1   %Yksilön siirtäminen toiseen populaatioon.
            inds = 1:ninds;
            aputaulu = [inds' rand(ninds,1)];
            aputaulu = sortrows(aputaulu,2);
            inds = aputaulu(:,1)';
            
            muutosNyt = 0;
            for ind = inds
                i1 = PARTITION_IN(ind);
                [muutokset, diffInCounts] = laskeMuutokset(ind, rowsFromInd, ...
                    data, adjprior, priorTerm);
                
                if round==1, [maxMuutos, i2] = max(muutokset); end
               
                if (i1~=i2 & maxMuutos>1e-5)
                    % Tapahtui muutos
                    muutoksia = 1;
                    kivaluku = kivaluku+1;
                    updateGlobalVariables(ind, i2, rowsFromInd, diffInCounts,...
                        adjprior, priorTerm);
                    logml = logml+maxMuutos;
                       
                end
            end
       
        elseif round==2  %Populaation yhdistäminen toiseen.
            maxMuutos = 0;
            for pop = 1:npops
                [muutokset, diffInCounts] = laskeMuutokset2(pop, rowsFromInd, ...
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
                updateGlobalVariables2(i1,i2,rowsFromInd, diffInCountsBest, ...
                    adjprior, priorTerm);
                logml = logml + maxMuutos;  
            end
                        
         
        elseif round==3 | round==4 %Populaation jakaminen osiin.
            maxMuutos = 0;
            ninds = size(data,1)/rowsFromInd;
            for pop = 1:npops
                inds2 = find(PARTITION_IN==pop);
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
                    muutokset = laskeMuutokset3(T2, inds2, rowsFromInd, data, ...
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
                rows = computeRows(rowsFromInd, muuttuvat, length(muuttuvat));
                diffInCounts = computeDiffInCounts(rows, size(COUNTS_IN,1), ...
                    size(COUNTS_IN,2), data);
                i1 = PARTITION_IN(muuttuvat(1));
                updateGlobalVariables3(muuttuvat, rowsFromInd, diffInCounts, ...
                    adjprior, priorTerm, i2);
                logml = logml + maxMuutos;
               
            end
    
        elseif round == 5 | round == 6
            pop=0;
            muutettu = 0;
            poplogml = POP_LOGML_IN;
            partition = PARTITION_IN;
            counts = COUNTS_IN;
            sumcounts = SUMCOUNTS_IN;
           
            while (pop < npops & muutettu == 0)
                pop = pop+1;
                totalMuutos = 0;
                inds = find(PARTITION_IN==pop);
                if round == 5
                    aputaulu = [inds rand(length(inds),1)];
                    aputaulu = sortrows(aputaulu,2);
                    inds = aputaulu(:,1)';
                elseif round == 6
                    inds = returnInOrder(inds, pop, rowsFromInd, data, adjprior, priorTerm);
                end
                
                i=0;
                           
                while (length(inds)>0 & i<length(inds))
                    i = i+1;
                    ind = inds(i);
                    [muutokset, diffInCounts] = laskeMuutokset(ind, rowsFromInd, ...
                        data, adjprior, priorTerm);
                    muutokset(pop) = -1e50;   % Varmasti ei suurin!!!
                    [maxMuutos, i2] = max(muutokset);
                    updateGlobalVariables(ind, i2, rowsFromInd, diffInCounts,...
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
                    muutettu=1;
                    muutoksia = 1;  % Ulompi kirjanpito.
                else
                    % Missään vaiheessa tila ei parantunut.
                    % Perutaan kaikki muutokset.
                    PARTITION_IN = partition;
                    SUMCOUNTS_IN = sumcounts;
                    POP_LOGML_IN = poplogml;
                    COUNTS_IN = counts;
                    logml = logml - totalMuutos;
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
            roundTypes = [2];
        elseif vaihe==3
            roundTypes=[5];
        elseif vaihe==4
            roundTypes=[4 3 1];
        elseif vaihe
            roundTypes=[6 2 3 4 1];
        end
    end
end

partition = PARTITION_IN;
counts = COUNTS_IN;
sumcounts = SUMCOUNTS_IN;



%-------------------------------------------------------------------------------------

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

initializeGammaln(ninds, rowsFromInd, max(noalle));

logml = computeLogml(counts, sumcounts, noalle, data, rowsFromInd);


%-----------------------------------------------------------------------


function logml=computeLogml(counts, sumcounts, noalle, data, rowsFromInd)
nloci = size(counts,2);
npops = size(counts,3);
adjnoalle = zeros(max(noalle),nloci);
for j=1:nloci
    adjnoalle(1:noalle(j),j)=noalle(j);
    if (noalle(j)<max(noalle))
        adjnoalle(noalle(j)+1:end,j)=1;
	end       
end

%logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
%    - npops*sum(sum(gammaln(adjprior))) - ...
%    sum(sum(gammaln(1+sumcounts)));
%logml = logml2;

global GAMMA_LN;
rowsInG = size(data,1)+rowsFromInd;

logml = sum(sum(sum(GAMMA_LN(counts+1 + repmat(rowsInG*(adjnoalle-1),[1 1 npops]))))) ...
    - npops*sum(sum(GAMMA_LN(1, adjnoalle))) ...
    -sum(sum(GAMMA_LN(sumcounts+1,1)));


%--------------------------------------------------------------------------


function initializeGammaln(ninds, rowsFromInd, maxAlleles)
%Alustaa GAMMALN muuttujan s.e. GAMMALN(i,j)=gammaln((i-1) + 1/j)
global GAMMA_LN;
GAMMA_LN = zeros((1+ninds)*rowsFromInd, maxAlleles);
for i=1:(ninds+1)*rowsFromInd
    for j=1:maxAlleles
        GAMMA_LN(i,j)=gammaln((i-1) + 1/j);
    end
end

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------

function rows = computeRows(rowsFromInd, inds, ninds)
% On annettu yksilöt inds. Funktio palauttaa vektorin, joka
% sisältää niiden rivien numerot, jotka sisältävät yksilöiden
% dataa.

rows = inds(:, ones(1,rowsFromInd));
rows = rows*rowsFromInd;
miinus = repmat(rowsFromInd-1 : -1 : 0, [ninds 1]);
rows = rows - miinus;
rows = reshape(rows', [1,rowsFromInd*ninds]);

%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, rowsFromInd, diffInCounts, ...
    adjprior, priorTerm)
% Suorittaa globaalien muuttujien muutokset, kun yksilö ind
% on siirretään koriin i2.

global PARTITION_IN; 
global COUNTS_IN; 
global SUMCOUNTS_IN;
global POP_LOGML_IN;

i1 = PARTITION_IN(ind);
PARTITION_IN(ind)=i2;

COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1) - diffInCounts;
COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2) + diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:) - sum(diffInCounts);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:) + sum(diffInCounts);

POP_LOGML_IN([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%---------------------------------------------------------------------------------


function updateGlobalVariables2( ...
    i1, i2, rowsFromInd, diffInCounts, adjprior, priorTerm);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilöt siirretään koriin i2.

global PARTITION_IN; 
global COUNTS_IN; 
global SUMCOUNTS_IN;
global POP_LOGML_IN;

inds = find(PARTITION_IN==i1);
PARTITION_IN(inds) = i2;

COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1) - diffInCounts;
COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2) + diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:) - sum(diffInCounts);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:) + sum(diffInCounts);

POP_LOGML_IN(i1) = 0;
POP_LOGML_IN(i2) = computePopulationLogml(i2, adjprior, priorTerm);


%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, rowsFromInd, diffInCounts, ...
    adjprior, priorTerm, i2);
% Suorittaa globaalien muuttujien päivitykset, kun yksilöt 'muuttuvat'
% siirretään koriin i2. Ennen siirtoa yksilöiden on kuuluttava samaan
% koriin.

global PARTITION_IN;
global COUNTS_IN;
global SUMCOUNTS_IN;
global POP_LOGML_IN;

i1 = PARTITION_IN(muuttuvat(1));
PARTITION_IN(muuttuvat) = i2;

COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1) - diffInCounts;
COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2) + diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:) - sum(diffInCounts);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:) + sum(diffInCounts);

POP_LOGML_IN([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%----------------------------------------------------------------------


function inds = returnInOrder(inds, pop, rowsFromInd, data, adjprior, priorTerm)
% Palauttaa yksilöt järjestyksessä siten, että ensimmäisenä on
% se, jonka poistaminen populaatiosta pop nostaisi logml:n
% arvoa eniten.

global COUNTS_IN;      global SUMCOUNTS_IN;
ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind = inds(i);
    rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
    diffInCounts = computeDiffInCounts(rows, size(COUNTS_IN,1), size(COUNTS_IN,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS_IN(:,:,pop) = COUNTS_IN(:,:,pop)-diffInCounts;
    SUMCOUNTS_IN(pop,:) = SUMCOUNTS_IN(pop,:)-diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
    COUNTS_IN(:,:,pop) = COUNTS_IN(:,:,pop)+diffInCounts;
    SUMCOUNTS_IN(pop,:) = SUMCOUNTS_IN(pop,:)+diffInSumCounts;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, rowsFromInd, data, adjprior, priorTerm)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikä olisi
% muutos logml:ssä, mikäli yksilö ind siirretään koriin i.
% diffInCounts on poistettava COUNTS_IN:in siivusta i1 ja lisättävä
% COUNTS_IN:in siivuun i2, mikäli muutos toteutetaan.

global COUNTS_IN;      global SUMCOUNTS_IN;
global PARTITION_IN;   global POP_LOGML_IN;
npops = size(COUNTS_IN,3);
muutokset = zeros(npops,1);

i1 = PARTITION_IN(ind);
i1_logml = POP_LOGML_IN(i1);

rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
diffInCounts = computeDiffInCounts(rows, size(COUNTS_IN,1), size(COUNTS_IN,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)-diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)+diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML_IN(i2);

COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2( ...
    i1, rowsFromInd, data, adjprior, priorTerm);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikä olisi
% muutos logml:ssä, mikäli korin i1 kaikki yksilöt siirretään
% koriin i. 

global COUNTS_IN;      global SUMCOUNTS_IN;
global PARTITION_IN;   global POP_LOGML_IN;
npops = size(COUNTS_IN,3);
muutokset = zeros(npops,1);

i1_logml = POP_LOGML_IN(i1);

inds = find(PARTITION_IN==i1);
ninds = length(inds);

if ninds==0
    diffInCounts = zeros(size(COUNTS_IN,1), size(COUNTS_IN,2));
    return;
end

rows = computeRows(rowsFromInd, inds, ninds);

diffInCounts = computeDiffInCounts(rows, size(COUNTS_IN,1), size(COUNTS_IN,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)-diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)+diffInCounts;
SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML_IN(i2);

COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;



%------------------------------------------------------------------------------------


function muutokset = laskeMuutokset3(T2, inds2, rowsFromInd, ...
    data, adjprior, priorTerm, i1)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mikä olisi muutos logml:ssä, jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirretään koriin j.

global COUNTS_IN;      global SUMCOUNTS_IN;
global PARTITION_IN;   global POP_LOGML_IN;
npops = size(COUNTS_IN,3);
npops2 = length(unique(T2));
muutokset = zeros(npops2, npops);

i1_logml = POP_LOGML_IN(i1);

for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
        rows = computeRows(rowsFromInd, inds, ninds);
        diffInCounts = computeDiffInCounts(rows, size(COUNTS_IN,1), size(COUNTS_IN,2), data);
        diffInSumCounts = sum(diffInCounts);

        COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)-diffInCounts;
        SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)-diffInSumCounts;
        new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
        COUNTS_IN(:,:,i1) = COUNTS_IN(:,:,i1)+diffInCounts;
        SUMCOUNTS_IN(i1,:) = SUMCOUNTS_IN(i1,:)+diffInSumCounts;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML_IN(i2)';
    
        COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
        new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
        COUNTS_IN(:,:,i2) = COUNTS_IN(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS_IN(i2,:) = SUMCOUNTS_IN(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

        muutokset(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end    
end


%------------------------------------------------------------------------------------

function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukumäärät (vastaavasti kuin COUNTS_IN:issa), jotka ovat data:n 
% riveillä rows.

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
% logml:t koreille, jotka on määritelty pops-muuttujalla.

global COUNTS_IN;
global SUMCOUNTS_IN;
x = size(COUNTS_IN,1);
y = size(COUNTS_IN,2);
z = length(pops);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS_IN(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS_IN(pops,:)),2) - priorTerm;


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