function c = preprocAln(alnMat)
% This function preprocess the alignment matrix to cliques and separators
% Lu Cheng, 24.05.2011

[nSeq nLoci] = size(alnMat);

alnCell = mat2cell(alnMat,nSeq,ones(1,nLoci));

arrUniqBase = cellfun(@unique,alnCell,'UniformOutput',false);  % unique base at each loci
arrUniqBaseNum = cellfun(@length,arrUniqBase);

arrCqNum = arrUniqBaseNum(1:end-2).*arrUniqBaseNum(2:end-1).*arrUniqBaseNum(3:end);
arrSpNum = arrUniqBaseNum(2:end-2).*arrUniqBaseNum(3:end-1);

nMaxCqCodes = max(arrCqNum);
nMaxSpCodes = max(arrSpNum);

cqCodes = cellfun(@myProd,arrUniqBase(1:end-2),arrUniqBase(2:end-1),arrUniqBase(3:end), ...
    'UniformOutput',false);
spCodes = cellfun(@myProd,arrUniqBase(2:end-2),arrUniqBase(3:end-1), ...
    'UniformOutput',false);

cqData = zeros(nSeq,length(cqCodes));
spData = zeros(nSeq,length(spCodes));

cqCounts = zeros(nMaxCqCodes,length(cqCodes));
spCounts = zeros(nMaxSpCodes,length(spCodes));

cqPrior = ones(nMaxCqCodes,length(cqCodes));
spPrior = ones(nMaxSpCodes,length(spCodes));

for i=1:nLoci-2
    
    nCodeTmp =  size(cqCodes{i},1);
    for j=1:nCodeTmp
        tmpInds = ismember(alnMat(:,i:i+2),cqCodes{i}(j,:),'rows');
        cqData(tmpInds,i) = j;
        cqCounts(j,i) = sum(tmpInds);
    end
    
    cqPrior(1:nCodeTmp,i) = 1/nCodeTmp;
    
    if i==1
        continue;
    end
    
    k=i-1;
    nCodeTmp = size(spCodes{k},1);
    for j=1:nCodeTmp
        tmpInds = ismember(alnMat(:,i:i+1),spCodes{k}(j,:),'rows');
        spData(tmpInds,k) = j;
        spCounts(j,k) = sum(tmpInds);
    end
    
    spPrior(1:nCodeTmp,k) = 1/nCodeTmp;
end

c.nSeq = nSeq;
% c.alnMat = alnMat;

c.arrUniqBase = arrUniqBase;
c.arrUniqBaseNum = arrUniqBaseNum;

c.nMaxCqCodes = nMaxCqCodes;
c.nMaxSpCodes = nMaxSpCodes;

c.cqCodes = cqCodes;
c.spCodes = spCodes;

c.cqData = cqData;
c.spData = spData;

c.cqCounts = cqCounts;
c.spCounts = spCounts;

c.cqPrior = cqPrior;
c.spPrior = spPrior;


function y = myProd(varargin)
% calculate the cartesian product for the input
% Lu Cheng, 24.05.2011

if nargin==2
    set1 = varargin{1};
    set2 = varargin{2};
    [t1 t2] = meshgrid(set1,set2);
    y = [t1(:) t2(:)];
elseif nargin==3
    set1 = varargin{1};
    set2 = varargin{2};
    set3 = varargin{3};
    [t1 t2 t3] = meshgrid(set1,set2,set3);
    y = [t1(:) t2(:) t3(:)];
else
    y = [];
end