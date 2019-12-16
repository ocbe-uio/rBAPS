function logml = computeTotalLogml
% compute the log marginal likelihood of the data
% Lu Cheng, 15.12.2012

global LOGML_TABLE;
global LOC_CQ_COUNTS;
global LOC_SP_COUNTS;


notEmpty = any(LOC_CQ_COUNTS,1);
npops = length(find(notEmpty == 1));

% the following codes added by Lu Cheng, 15.12.2012
% some lines might all be zero if some sequence is deleted
tmpIndsCq = find(any(LOC_CQ_COUNTS,2));
tmpIndsSp = find(any(LOC_SP_COUNTS,2));

locCqCounts = LOC_CQ_COUNTS(tmpIndsCq,notEmpty);
locSpCounts = LOC_SP_COUNTS(tmpIndsSp,notEmpty);

sumcliq=sum(locCqCounts, 2);
sumsep=sum(locSpCounts, 2);

ncliq = length(tmpIndsCq);
nsep = length(tmpIndsSp);
cliqsizes = sum(locCqCounts, 2)';
sepsizes = sum(locSpCounts, 2)';
cliqsizes = min([cliqsizes; npops*ones(1,ncliq)])';
sepsizes = min([sepsizes; npops*ones(1,nsep)])';

klikkitn = sum(sum(gammaln(locCqCounts + repmat(1./cliqsizes, [1 npops])))) ... 
                    - sum(npops*(gammaln(1./cliqsizes))) ...
                    - sum(gammaln(sumcliq + 1));
                                    
septn = sum(sum(gammaln(locSpCounts + repmat(1./sepsizes, [1 npops])))) ...
                - sum(npops*(gammaln(1./sepsizes))) ...
                - sum(gammaln(sumsep + 1));
            
spatialPrior = (klikkitn - septn);


logml = sum(LOGML_TABLE) + spatialPrior;