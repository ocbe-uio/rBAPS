function [counts,noalle_est,prior,adjprior,rawalleles] = allfreqsnew2(rawdata, noalle_est)
% Filename: allfreqsnew2.m
% [counts,noalle,prior,adjprior,rawalleles] = allfreqsnew(rawdata)
%
% Description:
% rawdata has n rows (2 x #individuals) and n(l) first
% colums are the loci. The last column is the subpopindex
% prior is a created matrix of positive Dirichlet hyperparameters
% missing data is filtered out
% !!!NEW!!!zeros are accepted as allele codes and any negative numbers as missing data.

% Modified by: Jing Tang
SCALE = 1;
dime=size(rawdata);
noalle=zeros(dime(2)-1,1);
rawalleles=cell(1,dime(2)-1);
for i=1:dime(2)-1
  noalle(i)=length(unique(rawdata(:,i)));
end
for i=1:dime(2)-1
   if length(find(rawdata(:,i)<=0))>0
      noalle(i)=noalle(i)-1;
   end
end

% Fomulate the raw data such that the value i in a entry denotes the ith
% alleles.
for i=1:dime(2)-1
    rawalles=unique(rawdata(:,i));
    % rawalles = [1:noalle(i)]';
    if rawalles(1)<=0
        rawalles(1)=-999;
    end
    rawalleles{i} = rawalles;    %rawalleles!!!
    if rawalles(1)<0
        for j=2:noalle(i)+1
          %rawdata(find(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*(j-1);
           rawdata(logical(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*(j-1);
        end
    else
        for j=1:noalle(i)
            % rawdata(find(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*j;
        rawdata(logical(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*j;
        end
    end
end

% ALLOWED_MEMORY = 50; % in unit of megabyte.
% n1 = max(noalle_est);
% n2 = dime(2)-1;
% n3 = double(max(rawdata(:,dime(2))));
% ncells = n1*n2*n3;
% memory_used = ncells/(1024*1024); % using uint8 format.
% if memory_used < ALLOWED_MEMORY
%     counts=zeros(n1,n2,n3,'uint8');
% else
%     nbatches = ceil(memory_used/ALLOWED_MEMORY);
%     n3_in_batch = ceil(n3/nbatches);
%     counts = cell(nbatches,1);
%       for i=1:nbatches-1
%           % counts = cat(3,counts,uint16(zeros(n1,n2,n3_in_batch)));
%           counts{i} = zeros(n1,n2,n3_in_batch,'uint8');
%       end
%       % counts = cat(3, counts, uint16(zeros(n1,n2,n3-n3_in_batch*(nbatches-1))));
%       counts{i} = zeros(n1,n2,n3-n3_in_batch*(nbatches-1),'uint8');
% end



counts = zeros(max(noalle_est),dime(2)-1,max(rawdata(:,dime(2))),'uint8');
for i=1:dime(1)
    for j=1:dime(2)-1
        if rawdata(i,j)>0
            counts(rawdata(i,j),j,rawdata(i,dime(2)))=...
                counts(rawdata(i,j),j,rawdata(i,dime(2)))+1;
        end
    end
end

maxnoalle = max(noalle_est);
% prior = [];
prior=zeros(maxnoalle,dime(2)-1);
for i=1:dime(2)-1
    prior(:,i) = [SCALE*ones(noalle_est(i),1)/noalle_est(i);zeros(maxnoalle-noalle_est(i),1)];
end

adjprior=prior;
for i=1:dime(2)-1
    adjprior(:,i)=adjprior(:,i)+[zeros(noalle_est(i),1);ones(maxnoalle-noalle_est(i),1)];
end
