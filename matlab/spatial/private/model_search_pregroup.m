function [partition, logml, partitionSummary, logmldiff] = model_search_pregroup(c, pgPart, pgDist, roundTypes, nMaxPops)
% This function clusters DNA alignment using "codon" model in Corander and Tang's
% paper: Bayesian analysis of population structure based on linked
% molecular information (2007), Mathematical Biosciences
% c: preprocessed data for the sequence alignment
% pgPart: partition which assign sequences to pregroups
% pgDist: distances between the pregroups
%       (1,2)(1,3)(1,4)...(2,3)(2,4).....(3,4)...(n-1,n)
% roundTypes: array of operation types

% Lu Cheng
% 21.03.2012

interactive = false;

global PARTITION; 
global CQ_COUNTS;global SUM_CQ_COUNTS; 
global SP_COUNTS;global SUM_SP_COUNTS; 
global CQ_PRIOR; global SP_PRIOR;
global LOGML_TABLE;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;

global LOC_SP_COUNTS;
global LOC_CQ_COUNTS;

clearGlobalVars;

nINDS = c.nSeq;
nPOPS = nMaxPops;

% load pregroup information
nPregroup = length(unique(pgPart));
if nPregroup<nMaxPops
    error('#pregroup: %d, nMaxPops: %d. Number of pregroups should be higher than maximum number of population. \n',nPregroup,nMaxPops);
end

pregroups = cell(nPregroup,1);
pgSize = zeros(nPregroup,1);
for i=1:nPregroup
    pregroups{i} = find(pgPart==i);
    pgSize(i) = length(pregroups{i});
end

pgZ = linkage(pgDist(:)','complete');
initPart = cluster(pgZ,'maxclust',nPOPS);
partition = zeros(nINDS,1);
for i=1:nPregroup
    partition(pregroups{i}) = initPart(i);
end
clear i pgZ initPart

% PRIOR VALUES:
CQ_PRIOR = c.cqPrior;
SP_PRIOR = c.spPrior;

% Initialize PARTITION, **_COUNTS, SUM_**_COUNTS, alnMat
[sumCqCounts, cqCounts] = initialCounts2(partition, c.cqData, nPOPS, c.nMaxCqCodes);
[sumSpCounts, spCounts] = initialCounts2(partition, c.spData, nPOPS, c.nMaxSpCodes);

CQ_COUNTS = cqCounts; SUM_CQ_COUNTS = sumCqCounts;
SP_COUNTS = spCounts; SUM_SP_COUNTS = sumSpCounts;

PARTITION = partition;
[cliqcounts, sepcounts] = computeCounts(c.locCliques, c.locSeparators, nPOPS);
LOC_CQ_COUNTS = cliqcounts;
LOC_SP_COUNTS = sepcounts;

partitionSummary = -Inf*ones(30,2,nPOPS);  % Tiedot 30 parhaasta partitiosta (npops ja logml)
partitionSummary(:,1,:) = zeros(30,1,nPOPS);
worstLogml = -Inf*ones(1, nPOPS); worstIndex = ones(1, nPOPS);

clear partition cqCounts sumCqCounts spCounts sumSpCounts;

% Initialize LOGML_TABLE:
nINDS = c.nSeq;
LOGML_TABLE = zeros(nPOPS,1);
updateLogmlTable(1:nPOPS);

REMOVAL_DIFFERENCE = zeros(nINDS,1);
REMOVAL_DIFFERENCE(:,:) = nan;
ADDITION_DIFFERENCE = zeros(nINDS,nPOPS);
ADDITION_DIFFERENCE(:,:) = nan;
JOIN_DIFFERENCE = zeros(nPOPS, nPOPS);
JOIN_DIFFERENCE(:,:) = nan;

% ***********Doc:********************
% REMOVAL_DIFFERENCE(ind) tells the change in logml if ind is removed from
% its cluster. nan, if the cluster has changed, since the value was last
% calculated.
% 
% ADDITION_DIFFERENCE(ind, pop) tells the change in logml if ind is added
% to cluster pop. nan, if the cluster has changed since the value was last
% calculated. Always nan, if pop is ind's own cluster.
%
% JOIN_DIFFERENCE(pop1,pop2) = tells the change in logml if pop1 and pop2
% are combined. nan, if either cluster has changed since the value was last
% calculated.
% ***********Doc end*****************

logml = computeTotalLogml;

disp('The beginning:');
% disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logml)]);
disp(' ');


nnotEmptyPops = length(unique(PARTITION));
if logml>worstLogml(nnotEmptyPops);
    [partitionSummary(:,:,nnotEmptyPops), added] = addToSummary(logml, ...
        partitionSummary(:,:,nnotEmptyPops), worstIndex(nnotEmptyPops));
    if (added==1)
        [worstLogml(nnotEmptyPops), worstIndex(nnotEmptyPops)] = ...
            min(partitionSummary(:,2,nnotEmptyPops));  
    end
end

% START SEARCH OF THE BEST PARTITION:

vipu = zeros(1,14); 
if interactive
    roundTypes = input('Input steps: ');
    if ischar(roundTypes), roundTypes = str2num(roundTypes); end
end
ready = 0;


while ready ~= 1

%     disp(['Performing steps: ' num2str(roundTypes)]);

    for n = 1:length(roundTypes)
        round = roundTypes(n);
        moveCounter = 0;

        if  round==1 && vipu(1)==0  % move an individual to another population
            
            pgInds = getMoveInds(pgPart,pgDist,nPregroup); % get pregroup inds to be moved

            for pgind = pgInds(:)'
%                 inds = cell2mat(pregroups(pgInds));
                tmpInds = pregroups{pgind};
                tmpChanges = calcLogmlChanges(tmpInds, c.cqData, c.nMaxCqCodes, ...
                            c.spData, c.nMaxSpCodes,  c.locCliques, c.locSeparators, logml);

                [maxChange, maxIndex] = max(tmpChanges);
                if maxChange>1e-5
                    updateGlobalVariables(tmpInds, maxIndex, c.cqData, c.nMaxCqCodes, ...
                        c.spData, c.nMaxSpCodes,c.locCliques, c.locSeparators);
%                     fprintf('moving from %d to %d.\n',PARTITION(ind),maxIndex)
                    logml = computeTotalLogml();
                    moveCounter = moveCounter+length(pgInds);
                    vipu = zeros(1,14);
                    
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
            if moveCounter==0, vipu(1)=1; end
            disp(['Step 1: ' num2str(moveCounter) ' pregroups were moved.']);
            
        elseif round==2 && vipu(2)==0  % join two populations

            update_join_difference(c.cqData, c.nMaxCqCodes, ...
                c.spData, c.nMaxSpCodes, c.locCliques, c.locSeparators, logml);
            [maxChange, aux] = max(JOIN_DIFFERENCE(:));
            [i1, i2] = ind2sub([nPOPS,nPOPS],aux);

            if maxChange>1e-5
                tmpInds = find(PARTITION==i1);
                updateGlobalVariables(tmpInds, i2, c.cqData, c.nMaxCqCodes, ...
                        c.spData, c.nMaxSpCodes, c.locCliques, c.locSeparators);
                logml = computeTotalLogml;

                disp(['Step 2: Clusters ' num2str(i1) ' and ' num2str(i2) ' combined.']);
                vipu = zeros(1,14);
                
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
                disp('Step 2: no changes.');
                vipu(2)=1;
            end
        elseif ismember(round, 3:4) && vipu(round)==0  % Split a population, and move one subpopulation to another population

            pops = randperm(nPOPS);

            splitFlags = zeros(nPOPS,1);
            for pop = pops(:)'

                maxChange = 0;
                indsToBeMoved = [];

                inds2 = find(PARTITION==pop);
                pgInds2 = unique(pgPart(inds2));
                nPgInds2 = length(unique(pgPart(inds2)));
                if nPgInds2>4

                    if round==3
                        dist3 = getDistance(pgInds2,pgDist,nPregroup);
                        npops2 = min(20, floor(nPgInds2 / 5));
                    elseif round==4
                        dist3 = getDistance(pgInds2,pgDist,nPregroup);
                        npops2 = 2;
                    end

                    Z3 = linkage(dist3(:)','complete');
                    T3 = cluster(Z3, 'maxclust', npops2);

                    for i = 1:npops2
                        indsX = pgInds2(T3==i);
                        indsX = cell2mat(pregroups(indsX));
                        tmpChanges = calcLogmlChanges(indsX, c.cqData, c.nMaxCqCodes, ...
                            c.spData, c.nMaxSpCodes, c.locCliques, c.locSeparators, logml);
                        [tmpMaxChange, tmpMaxPop] = max(tmpChanges);
                        if tmpMaxChange>maxChange
                            maxChange = tmpMaxChange;
                            % i1 = pop;
                            i2 = tmpMaxPop;
                            indsToBeMoved = indsX;
                        end
                    end
                    if maxChange>1e-5
                        updateGlobalVariables(indsToBeMoved, i2, c.cqData, c.nMaxCqCodes, ...
                            c.spData, c.nMaxSpCodes, c.locCliques, c.locSeparators);
                        logml = computeTotalLogml;
                        splitFlags(pop)=1;
                        
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
            end
            if any(splitFlags)
                disp(['Step ' num2str(round) ': ' num2str(sum(splitFlags)) ' populations were split.']);
                vipu = zeros(1,14);
            else
                disp(['Step ' num2str(round) ': no changes.']);
                vipu(round)=1;
            end
        end        
    end

    if interactive
        roundTypes = input('Input extra steps: ');
        if ischar(roundTypes), roundTypes = str2num(roundTypes); end
    else
        roundTypes = [];
    end

    if isempty(roundTypes)
        ready = 1;
    end
end

% disp(' ');
disp('BEST PARTITION: ');
% disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml): ' num2str(logml)]);
disp(' ');

nPOPS = rmEmptyPopulation(c.locCliques, c.locSeparators);

logmldiff = zeros(nPregroup,nPOPS);  % the change of logml if pregroup i is moved to group j
for i=1:nPregroup
    tmpInds = pregroups{i};
    tmpChanges = calcLogmlChanges(tmpInds, c.cqData, c.nMaxCqCodes, ...
                            c.spData, c.nMaxSpCodes,  c.locCliques, c.locSeparators, logml);
    logmldiff(i,:) = tmpChanges';                   
end
logmldiff(isnan(logmldiff))=0;

partition = zeros(nPregroup,1);
for i=1:nPregroup
    partition(i)=unique(PARTITION(pgPart==i));
end

%----------------------------------------------------------------------------


function [dist2, dind1, dind2] = getDistance(inds2, origDist, ninds)
% pick out the distrances between samples in "inds2" from "origDist"
% origDist specifies the distances of (1,2),(1,3),(1,4)......(ninds-1,ninds)
% Lu Cheng, 22.06.2011

if ~issorted(inds2)
    error('inds2 is not in ascending order!');
end

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
irow = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(irow, 1) = inds2(i);
        apu(irow, 2) = inds2(j);
        irow = irow+1;
    end
end

dind1 = apu(:,1);
dind2 = apu(:,2);

apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = origDist(apu);

%---------------------------------------------------------------


function inds = getMoveInds(pgPart, pgDist, nPregroup)
% get pregroup indexs to be moved to another cluster
% we always take the 35% pregroups of each cluster which are most distant
% to each other
% Lu Cheng, 22.06.2011

global PARTITION;

pops = unique(PARTITION);
inds = [];

for tmpPop = pops(:)'
    tmpInds = unique(pgPart(PARTITION==tmpPop));

    if(length(tmpInds)<20)
        inds = [inds tmpInds(:)']; %#ok<AGROW>
        continue;
    end
    
    [tmpDist, dind1, dind2] = getDistance(tmpInds,pgDist,nPregroup);
    tmpVal = quantile(tmpDist,0.65);
    tmpInds2 = find(tmpDist>tmpVal);
    tmpInds3 = union(unique(dind1(tmpInds2)), unique(dind2(tmpInds2)));
    inds = [inds tmpInds3(:)']; %#ok<AGROW>
end


% ------------------------------------------------------------

function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedet‰‰n, ett?annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss?ei viel?ole
% annettua logml arvoa, niin lis‰t‰‰n worstIndex:in kohtaan uusi logml ja
% nykyist?partitiota vastaava nclusters:in arvo. Muutoin ei tehd?mit‰‰n.
global PARTITION;

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lˆydetty partitio ei ole viel?kirjattuna summaryyn.

    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end




