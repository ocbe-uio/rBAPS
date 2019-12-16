function updateLogmlTable(pops)
% Updates global variables LOGML_TABLE, npops*1 array, logml values for
% each population given in "pops"
% After the updates, the values are based on the current values of the
% global variables CQ_COUNTS, SUM_CQ_COUNTS, SP_COUNTS, SUM_SP_COUNTS
% Lu Cheng, 25.05.2011

global CQ_COUNTS; global SUM_CQ_COUNTS; 
global SP_COUNTS; global SUM_SP_COUNTS; 
global CQ_PRIOR; global SP_PRIOR; 

global LOGML_TABLE;

tmpN = length(pops);
tmpCqPrior = repmat(CQ_PRIOR,[1 1 tmpN]);
tmpSpPrior = repmat(SP_PRIOR,[1 1 tmpN]);

term1 = 0-gammaln(1+SUM_CQ_COUNTS(:,pops));
term2 = sum(gammaln(tmpCqPrior+CQ_COUNTS(:,:,pops))-gammaln(tmpCqPrior) , 1);
if length(pops) > 1
    term2 = squeeze(term2);
else
    term2 = term2';
end

term3 = 0-gammaln(1+SUM_SP_COUNTS(:,pops));
term4 = sum(gammaln(tmpSpPrior+SP_COUNTS(:,:,pops))-gammaln(tmpSpPrior) , 1);

if length(pops) > 1
    term4 = squeeze(term4);
else
    term4 = term4';
end

LOGML_TABLE(pops) = sum(term1+term2) - sum(term3+term4);

%----------------------------------------------------------------------