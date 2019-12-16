function update_join_difference(cqData, nCqCodes, spData, nSpCodes, locCliques, locSeparators, logml)
% update JOIN_DIFFERENCE
% Lu Cheng, 15.12.2012

global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS;
global PARTITION;
global JOIN_DIFFERENCE;

global LOC_CQ_COUNTS;
global LOC_SP_COUNTS;

npops = size(CQ_COUNTS,3);

for i1 = 1:npops-1
    indsToBeMoved = find(PARTITION==i1);
    if isempty(indsToBeMoved)
        % Cluster i1 is empty
        JOIN_DIFFERENCE(i1,(i1+1):npops) = 0;
        JOIN_DIFFERENCE((i1+1):npops,i1) = 0;
    else
        [diffCqCounts diffCqSumCounts] = computeDiffInCounts(indsToBeMoved, cqData, nCqCodes);
        [diffSpCounts diffSpSumCounts] = computeDiffInCounts(indsToBeMoved, spData, nSpCodes);
        diffLocCqCounts = computeDiffInCliqCounts(locCliques, indsToBeMoved);
        diffLocSpCounts = computeDiffInCliqCounts(locSeparators, indsToBeMoved);
        
        unknown_pops = find(isnan(JOIN_DIFFERENCE(i1,(i1+1):end)));
        unknown_pops = unknown_pops+i1;
        
        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) - diffCqCounts;
        SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) - diffCqSumCounts;
        
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) - diffSpCounts;
        SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) - diffSpSumCounts;
        
        LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) - diffLocCqCounts;
        LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) - diffLocSpCounts;
        
%         PARTITION(indsToBeMoved) = -1;
        updateLogmlTable(i1);
        
        for i2 = unknown_pops
            CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) + diffCqCounts;
            SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) + diffCqSumCounts;
        
            SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) + diffSpCounts;
            SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) + diffSpSumCounts;
            
            LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) + diffLocCqCounts;
            LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) + diffLocSpCounts;
        
%             PARTITION(indsToBeMoved) = i2;
            updateLogmlTable(i2);
            logml_new = computeTotalLogml();
            JOIN_DIFFERENCE(i1,i2) = logml_new-logml;
            JOIN_DIFFERENCE(i2,i1) = logml_new-logml;
            
            CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2) - diffCqCounts;
            SUM_CQ_COUNTS(:,i2) = SUM_CQ_COUNTS(:,i2) - diffCqSumCounts;
        
            SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2) - diffSpCounts;
            SUM_SP_COUNTS(:,i2) = SUM_SP_COUNTS(:,i2) - diffSpSumCounts;
            
            LOC_CQ_COUNTS(:,i2) = LOC_CQ_COUNTS(:,i2) - diffLocCqCounts;
            LOC_SP_COUNTS(:,i2) = LOC_SP_COUNTS(:,i2) - diffLocSpCounts;
                        
%             PARTITION(indsToBeMoved) = -1;
            updateLogmlTable(i2);
        end
        
        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1) + diffCqCounts;
        SUM_CQ_COUNTS(:,i1) = SUM_CQ_COUNTS(:,i1) + diffCqSumCounts;
        
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1) + diffSpCounts;
        SUM_SP_COUNTS(:,i1) = SUM_SP_COUNTS(:,i1) + diffSpSumCounts;
        
        LOC_CQ_COUNTS(:,i1) = LOC_CQ_COUNTS(:,i1) + diffLocCqCounts;
        LOC_SP_COUNTS(:,i1) = LOC_SP_COUNTS(:,i1) + diffLocSpCounts;
        
%         PARTITION(indsToBeMoved) = i1;
        updateLogmlTable(i1);
    end
end