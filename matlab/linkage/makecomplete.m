function [data, totalmissing] = makecomplete(data)
%MAKECOMPLETE ESTIMATING missing alleles in a loci. 
%   Input: 
%     data: a baps data format, e.g. last column is index. Missing values
%     are denoted as any non-positive integers, commonly as 0 or -999.


[nrows, ncols] = size(data);
id = all(data>0);
if all(id) % no missing values
    return 
else
    fprintf(1,'Estimating the missing values...');
    
    % remove the totally missing loci
    missingloci = find(id==0);
    totalmissing = find(all(data(:,missingloci)<=0));
    if ~isempty(totalmissing) 
        disp('Totally missing loci were found.');
        data(:,missingloci(totalmissing)) = ones(nrows,length(totalmissing));
    end
%     goodloci = setdiff([1:ncols],missingloci(totalmissing));
%     data = data(:,goodloci);
%     missingloci = find(all(data>0)==0);
    
    data_in = data(:,[missingloci end]);
    [counts,noalle,prior,adjprior,rawalleles,data_in] = allfreqsnew3(data_in);
    data(:,missingloci) = data_in(:,[1:end-1]);
    
    s_counts = sum(counts,3);
%     m_counts = size(counts,3)*ones(1, size(counts,2)) - sum(s_counts); % counts of missing values
    m_counts = nrows*ones(1, size(counts,2)) - sum(s_counts); % counts of missing values
    m = repmat(m_counts,size(s_counts,1),1);
    a = s_counts + prior;
    a = a./repmat(sum(a),size(a,1),1);
    c = round(m.*a);
    
    % Be sure that the sum matches
    reassign = find(sum(c)~=m_counts);
    for i = reassign
        dif = sum(c(:,i)) - m_counts(:,i);
        if dif > 0
            remove = find(c(:,i)==max(c(:,i)));
            c(remove(1),i)=c(remove(1),i) - dif; % decrease by the diff
        else
            add = find(c(:,i)==min(c(:,i)));
            c(add(1),i)=c(add(1),i) - dif; % increase by the diff
        end
    end
    try
    % Fill in the missing values    
    for i = 1:length(missingloci)
        missingrows = find(data(:,missingloci(i))<=0);
        p = randperm(m_counts(i));
        cum = [0 cumsum(c(:,i))'];
        for j = 1:size(c,1)
            data(missingrows(p([cum(j)+1:cum(j+1)])),missingloci(i)) = j;
        end
    end
    catch
        % disp('*** ERROR: in completing the missing data.');
        data = [];
        totalmissing = [];
        return
    end
    fprintf(1,'Finished.\n');
end



