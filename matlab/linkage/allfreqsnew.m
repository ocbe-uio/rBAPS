function [counts,noalle,prior,adjprior,rawalleles,rawdata] = allfreqsnew(rawdata)
% Filename: allfreqsnew.m
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
   if length(find(rawdata(:,i)<=0))>0 % filtering out the missing data
      noalle(i)=noalle(i)-1;
   end
end
for i=1:dime(2)-1
    rawalles=unique(rawdata(:,i));
    if rawalles(1)<=0
        rawalles(1)=-999;
    end
    rawalleles{i} = rawalles;    % rawalleles!!!
    if rawalles(1)<=0
    for j=2:noalle(i)+1
        rawdata(find(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*(j-1);
    end
    else
    for j=1:noalle(i)
        rawdata(find(rawdata(:,i)==rawalles(j)),i)=ones(length(find(rawdata(:,i)==rawalles(j))),1)*j;
    end
end
end

        
counts=zeros(max(noalle),dime(2)-1,max(rawdata(:,dime(2))));

for i=1:dime(1)
   for j=1:dime(2)-1
      if rawdata(i,j)>0
         counts(rawdata(i,j),j,rawdata(i,dime(2)))=...
            counts(rawdata(i,j),j,rawdata(i,dime(2)))+1;
      end
   end
end

prior=[];
for i=1:dime(2)-1
   prior=[prior [SCALE*ones(noalle(i),1)/noalle(i);zeros(max(noalle)-noalle(i),1)]];
end

adjprior=prior;
for i=1:dime(2)-1
   adjprior(:,i)=adjprior(:,i)+[zeros(noalle(i),1);ones(max(noalle)-noalle(i),1)];
end
