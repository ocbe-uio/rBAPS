function  compare(varargin)
% COMPARE compares the results from multiple runs.
%  input: is a group of result .mat files on the same data.

% Example: compare('e:\data\result1.mat','e:\data\result2.mat',...)
% or call it from the BAPS menu.

if nargin == 1
    error('number of input arguments must be >=2');
end

if nargin == 0
    out = uipickfiles('FilterSpec','*.mat',...
               'Prompt','Select mixture results: be sure that the underlying data and models are consistent.');
    if isnumeric(out)
        return
    end
    nfiles = length(out);
    filesin = out;

else
    nfiles = nargin;
    filesin = varargin;
end

display('---------------------------------------------------');
fprintf(1,'Comparing results ...\n');
if nfiles == 1 
    disp('*** ERROR: Too few files.');
    return
end
for i = 1:nfiles
    struct_array = load(filesin{i});
    if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        clear struct_array;
        if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd') || strcmp(c.mixtureType,'admix')
            fprintf(1,'*** ERROR: Incorrect mixture result in file %d\n',i );
            return
        end
    elseif isfield(struct_array,'PARTITION')  %Mideva versio
        c = struct_array;
        if ~isfield(c,'rowsFromInd')
            fprintf(1,'*** ERROR: Incorrect mixture result in file %d\n',i );
            return
        end
    else
        fprintf(1,'*** ERROR: Incorrect mixture result in file %d\n',i );
        return;
    end
    
    try
    partitionMat(i,:) = sort_partition(c.PARTITION);
    catch
        error('*** ERROR: inconsistent results.');
    end
    mixtureType{i} = c.mixtureType;
    logml(i) = c.logml;
    clear c;
end

len_mixture = length(mixtureType{1});
for i = 2:nfiles
    if len_mixture ~= length(mixtureType{i});
        error('*** ERROR: inconsistent mixture types.');
    end
end


% Find the best partition
best = logical(logml == max(logml));
[uniquepartition, ind1, ind2] = unique(partitionMat(best,:), 'rows');
fprintf(1,'Best partition was found at ''%s''\n',filesin{best});






    

