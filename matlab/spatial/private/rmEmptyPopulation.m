function [npops notEmpty] = rmEmptyPopulation(locCliques,locSeparators)
% remove empty populations from CQ_COUNTS and SUM_CQ_COUNTS, SP_COUNTS,
% SUM_SP_COUNTS
% update PARTITION
% Lu Cheng, 15.12.2012

global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS; 
global PARTITION;

global LOGML_TABLE;
global ADDITION_DIFFERENCE;
global JOIN_DIFFERENCE;

global LOC_CQ_COUNTS;
global LOC_SP_COUNTS;

notEmpty = find(any(SUM_CQ_COUNTS,1) & any(SUM_SP_COUNTS,1));

CQ_COUNTS = CQ_COUNTS(:,:,notEmpty);
SP_COUNTS = SP_COUNTS(:,:,notEmpty);

SUM_CQ_COUNTS = SUM_CQ_COUNTS(:,notEmpty);
SUM_SP_COUNTS = SUM_SP_COUNTS(:,notEmpty);

LOGML_TABLE = LOGML_TABLE(notEmpty);
ADDITION_DIFFERENCE = ADDITION_DIFFERENCE(:,notEmpty);
JOIN_DIFFERENCE = JOIN_DIFFERENCE(notEmpty,notEmpty);

for i=1:length(notEmpty)
    apu = (PARTITION==notEmpty(i));
    PARTITION(apu)=i;
end

npops = length(notEmpty);

[cliqcounts, sepcounts] = computeCounts(locCliques, locSeparators, npops);

LOC_CQ_COUNTS = cliqcounts;
LOC_SP_COUNTS = sepcounts;

