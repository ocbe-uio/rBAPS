function changes = calcLogmlChanges(inds, cqData, nCqCodes, spData, nSpCodes, locCliques, locSeparators, logml)
% compute the logml change if the given inds are moved to another cluster
% the input inds are supposed to come from the same cluster
% changes is a npops*1 vector
% Lu Cheng, 15.12.2012

global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS;
global PARTITION;

global LOC_CQ_COUNTS;
global LOC_SP_COUNTS;

npops = size(CQ_COUNTS,3);
changes = zeros(npops,1);
indsToBeMoved = inds;

if isempty(indsToBeMoved), return, end

i1 = PARTITION(indsToBeMoved(1));
[diffCqCounts diffCqSumCounts]= computeDiffInCounts(indsToBeMoved, cqData, nCqCodes);
[diffSpCounts diffSpSumCounts]= computeDiffInCounts(indsToBeMoved, spData, nSpCodes);

diffLocCqCounts = computeDiffInCliqCounts(locCliques, indsToBeMoved);
diffLocSpCounts = computeDiffInCliqCounts(locSeparators, indsToBeMoved);

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) - diffCqCounts;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) - diffSpCounts;

SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) - diffCqSumCounts;
SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) - diffSpSumCounts;

LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) - diffLocCqCounts;
LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) - diffLocSpCounts;

% PARTITION(inds) = -1;
updateLogmlTable(i1);

for i2 = 1:npops
    if i2 ~= i1
        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) + diffCqCounts;
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) + diffSpCounts;
    
        SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) + diffCqSumCounts;
        SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) + diffSpSumCounts;
        
        LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) + diffLocCqCounts;
        LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) + diffLocSpCounts;
        
%         PARTITION(inds) = i2;
        updateLogmlTable(i2);
        logml_new = computeTotalLogml();
        changes(i2) = logml_new - logml;
        
        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) - diffCqCounts;
        SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) - diffCqSumCounts;
        
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) - diffSpCounts;
        SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) - diffSpSumCounts;
        
        LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) - diffLocCqCounts;
        LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) - diffLocSpCounts;
        
%         PARTITION(inds) = -1;
        updateLogmlTable(i2);
    end
end

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) + diffCqCounts;
SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) + diffCqSumCounts;
        
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) + diffSpCounts;
SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) + diffSpSumCounts;

LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) + diffLocCqCounts;
LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) + diffLocSpCounts;

% PARTITION(inds) = i1;
updateLogmlTable(i1);


%---------------------------------------------------------------------