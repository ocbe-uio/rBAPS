function update_difference_tables(ind, cqData, nCqLetter, ...
    spData, nSpLetter, locCliques, locSeparators,logml)
% update ADDITION_DIFFERENCE and REMOVAL_DIFFERENCE
% Lu Cheng, 15.12.2012

global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS; 
global PARTITION;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;

global LOC_CQ_COUNTS;
global LOC_SP_COUNTS;

rem_old = REMOVAL_DIFFERENCE;
add_old = ADDITION_DIFFERENCE;

[diffCqCounts diffCqSumCounts] = computeDiffInCounts(ind, cqData, nCqLetter); 
[diffSpCounts diffSpSumCounts] = computeDiffInCounts(ind, spData, nSpLetter); 
diffLocCqCounts = computeDiffInCliqCounts(locCliques, ind);
diffLocSpCounts = computeDiffInCliqCounts(locSeparators, ind);

i1 = PARTITION(ind);

if isnan(rem_old(ind))
    % Update removal difference for the individual:
    % note that we did NOT add the removed item to other clusters
    CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) - diffCqCounts;
    SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) - diffSpCounts;
    
    SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) - diffCqSumCounts;
    SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) - diffSpSumCounts;
    
    LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) - diffLocCqCounts;
    LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) - diffLocSpCounts;
    
%     PARTITION(ind) = -1;
    updateLogmlTable(i1);
    logml_new = computeTotalLogml();
    rem_old(ind) = logml_new-logml;
    
    CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) + diffCqCounts;
    SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) + diffSpCounts;
    
    SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) + diffCqSumCounts;
    SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) + diffSpSumCounts;
    
    LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) + diffLocCqCounts;
    LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) + diffLocSpCounts;

%     PARTITION(ind) = i1;
    updateLogmlTable(i1);
end

new_pops = isnan(add_old(ind,:));
new_pops(i1) = 0;   % Own cluster needs never be calculated.
new_pops = find(new_pops);

for i2 = new_pops(:)'
    % Update addition differences for the individual:
    % note that we did NOT remove the item
    CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) + diffCqCounts;
    SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) + diffSpCounts;
    
    SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) + diffCqSumCounts;
    SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) + diffSpSumCounts;
    
    LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) + diffLocCqCounts;
    LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) + diffLocSpCounts;
    
%     PARTITION(ind) = i2;
    updateLogmlTable(i2);
    logml_new = computeTotalLogml();
    add_old(ind,i2) = logml_new - logml;
    
    CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) - diffCqCounts;
    SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) - diffSpCounts;
    
    SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) - diffCqSumCounts;
    SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) - diffSpSumCounts;

    LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) - diffLocCqCounts;
    LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) - diffLocSpCounts;
    
%     PARTITION(ind) = i1;
    updateLogmlTable(i2);
end

REMOVAL_DIFFERENCE = rem_old;
ADDITION_DIFFERENCE = add_old;

%---------------------------------------------------------------------