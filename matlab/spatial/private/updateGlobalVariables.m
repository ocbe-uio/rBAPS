function updateGlobalVariables(inds, i2, cqData, nCqCodes, spData, nSpCodes, locCliques, locSeparators)
% this function moves the samples specified by "inds" to cluser i2
% then update all the global variables, "inds" are supposed to come from the
% same cluster
% Lu Cheng, 15.12.2012

global PARTITION; 
global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS; 
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;

global LOC_SP_COUNTS;
global LOC_CQ_COUNTS;

i1 = PARTITION(inds(1));
PARTITION(inds)=i2;

[diffCqCounts diffCqSumCounts]= computeDiffInCounts(inds, cqData, nCqCodes);
[diffSpCounts diffSpSumCounts]= computeDiffInCounts(inds, spData, nSpCodes);

diffLocCqCounts = computeDiffInCliqCounts(locCliques, inds);
diffLocSpCounts = computeDiffInCliqCounts(locSeparators, inds);

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) - diffCqCounts;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) - diffSpCounts;

SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) - diffCqSumCounts;
SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) - diffSpSumCounts;

LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) - diffLocCqCounts;
LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) - diffLocSpCounts;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) + diffCqCounts;
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) + diffSpCounts;
    
SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) + diffCqSumCounts;
SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) + diffSpSumCounts;

LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) + diffLocCqCounts;
LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) + diffLocSpCounts;

updateLogmlTable([i1 i2]);

REMOVAL_DIFFERENCE(PARTITION==i1) = nan;
REMOVAL_DIFFERENCE(PARTITION==i2) = nan;
ADDITION_DIFFERENCE(:,[i1 i2]) = nan;

JOIN_DIFFERENCE(:,i2) = nan;
JOIN_DIFFERENCE(i2,:) = nan;

if ~any(PARTITION==i1)
    % i1 became empty
    JOIN_DIFFERENCE(:,i1) = 0;
    JOIN_DIFFERENCE(i1,:) = 0;
    JOIN_DIFFERENCE(i1,i1) = nan;
else
    JOIN_DIFFERENCE(:,i1) = nan;
    JOIN_DIFFERENCE(i1,:) = nan;
end