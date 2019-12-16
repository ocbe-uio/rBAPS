function  compare_admix(varargin)
% COMPARE compares the results from multiple runs for admixture results
%  input: is a group of result .mat files on the same data.

% Example: compare('e:\data\result1.mat','e:\data\result2.mat',...)
% or call it from the BAPS menu.

if nargin == 1
    error('number of input arguments must be >=2');
end

if nargin == 0
    out = uipickfiles('FilterSpec','*.mat',...
            'Prompt','Select admixture results: be sure that the underlying data and parameters are consistent.');
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
minsize = zeros(nfiles,1);
iters = zeros(nfiles,1);
refInds = zeros(nfiles,1);
refIters = zeros(nfiles,1);
prop = cell(nfiles,1);

clusters = cell(nfiles,1);
adjprior = [];

if nfiles == 1 
    disp('*** ERROR: Too few files.');
    return
end

% read admixture files
for i = 1:nfiles
    struct_array = load(filesin{i});
    if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        clear struct_array;
        if ~isfield(c,'PARTITION') || ~isfield(c,'rowsFromInd') ...
                 || ~isfield(c,'proportionsIt')
            fprintf(1,'*** ERROR: Incorrect admixture result in file %d\n',i );
            return
        end
    elseif isfield(struct_array,'PARTITION')  %Mideva versio
        c = struct_array;
        if ~isfield(c,'rowsFromInd')
            fprintf(1,'*** ERROR: Incorrect admixture result in file %d\n',i );
            return
        end
    else
        fprintf(1,'*** ERROR: Incorrect admixture result in file %d\n',i );
        return;
    end
    
    prop{i} = c.proportionsIt;
    pvalue(:,i) = c.pvalue;
    clusters{i} = c.clusters;
    popnames = c.popnames;
    
    % parameters
    minsize(i) = c.minsize;
    iters(i) = c.iters;
    refInds(i) = c.refInds;
    refIters(i) = c.refIters;
    
    if i==1 
        adjprior = c.adjprior;
    else
        if ~isequal(adjprior,c.adjprior)
            disp('*** ERROR: incosistent admixture results.');
            return
        end
    end
    
    clear c;
end

if length(unique(minsize))~=1 || length(unique(iters))~=1 ...
        || length(unique(refInds))~=1 || length(unique(refIters))~=1
    disp('*** ERROR: inconsistent admixture parameters.');
    return
end

% now combine the results
prop_combine = prop{1};
[ninds npops] = size(prop_combine);
[pvalue_combine,index] = min(pvalue,[],2);
for i = 1:ninds
prop_combine(i,:) = prop{index(i)}(i,:);
end

% display the results
tulostaAdmixtureTiedot(prop_combine, pvalue_combine, minsize(1), iters(1)); 

viewPartition(prop_combine, popnames);

% save the results
talle = questdlg(['Do you want to save the combined admixture results?'], ...
    'Save results?','Yes','No','Yes');
if isequal(talle,'Yes')
    %waitALittle;
    [filename, pathname] = uiputfile('*.mat','Save results as');
   
    if (filename == 0) & (pathname == 0)
        % Cancel was pressed
        return
    else % copy 'baps4_output.baps' into the text file with the same name.
        if exist('baps4_output.baps','file')
            copyfile('baps4_output.baps',[pathname filename '.txt'])
            delete('baps4_output.baps')
        end
    end

    struct_array = load(filesin{1});
    c = struct_array.c;
    clear struct_array;
    c.proportionsIt = prop_combine;
    c.pvalue = pvalue_combine; % Added by Jing
    
    fprintf(1, 'Saving the results...');
%     save([pathname filename], 'c');
    save([pathname filename], 'c','-v7.3'); % added by Lu Cheng, 08.06.2012
    fprintf(1,'finished.\');

end






    

