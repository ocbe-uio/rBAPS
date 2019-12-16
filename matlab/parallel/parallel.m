function parallel(varargin)
% PARALLEL the main function of doing parallel classification/admixture.
%   input: order of input is arbitrary. The first option is the default.
%   'datafile' - the full path of the data.
%   'mixturetype', ['mix';'codon_mix';'linear_mix';'spatial';'ad_mix'];
%         - the classification/admixture model.
%   'initialk', a row vector of positive integers;
%         - the initial number of clusters.
%   'fixedk',['no';'yes'];
%         - whether the number of clusters is fixed during the computation.
%   'outputmat' - the full path of the output .mat file.
%   'datatype', ['numeric';'sequence';'matlab';'excel';'genepop']
%         - the data format;
%   'namefile' - the full path of the population name file.
%   'indexfile' - the full path of the index file.
%   'linkagemap' - the full path of the linkage map, needed only for the
%   unpreprocessed data under the linkage model.
%   'coordinatefile' - needed with the spatial model.
%   'groups', ['no';'yes'] - clustering of groups instead of individuals.
%   'groupname' - the full path of the group name file.

%   Examples:
%   - Linkage model:
%     parallel('datafile','e:\baps4\baps_source\data\bpseudomallei.xls',...
%            'mixturetype','codon_mix','initialk','[10 15]','fixedk','no',...
%            'outputmat','E:\test_link.mat','datatype','excel')
%   - Independent model:
%     parallel('datafile','e:\baps4\baps_source\data\baps_data.txt',...
%           'mixturetype','mix','initialk','[10:15]','fixedk','no',...
%           'outputmat','e:\test_ind.mat','datatype','numeric');
%   - Spatial model:
%     parallel('datafile','e:\baps4\baps_source\data\wolverines_spatial_preprocessed.mat',...
%              'mixturetype','spatial','initialk', '[10 11]', 'fixedk','no',...
%              'outputmat','e:\test_spatial.mat','datatype','matlab');
% 
%
%   - Admixture model:
%     parallel('datafile','e:\baps5\data\data1_mixture.mat', ...
%              'mixturetype','ad_mix', ...
%              'clusters','[1 3 5]',...
%              'iters','2',...
%              'refinds','3',...
%              'refiters','4',...
%              'outputmat','e:\baps5\data\data1_admixture_parallel.mat');

%   A group of result files can be later compared by using compare.m
%   function.
%-------------------------------------------------------------------------------
%- Set up options and default parameters
%-------------------------------------------------------------------------------
msgInvalidPair = '***ERROR: Bad value for argument: ''%s''';

% default options
options = struct('dataFile', '',...
    'dataType', 'numeric',...
    'mixtureType', 'mix',...
    'initialK', 1, ...
    'fixedK', 'no', ...
    'nameFile', '', ...
    'indexFile', '', ...
    'outputMat', '', ...
    'linkageMap','', ...
    'coordinateFile', '', ...
    'groups','no', ...
    'groupname', '', ...
    'clusters', '', ...
    'minSize', '', ...
    'iters', '', ...
    'refInds', '', ...
    'refIters', '' ...
    );

if nargin == 1 && isstruct(varargin{1})
    paramlist = [ fieldnames(varargin{1}) ...
        struct2cell(varargin{1}) ]';
    paramlist = { paramlist{:} };
else
    if mod(nargin,2)
        error('Invalid parameter/value pair arguments.');
    end
    paramlist = varargin;
end

optionsnames = lower(fieldnames(options));
for i=1:2:length(paramlist)
	pname = paramlist{i};
	pvalue = paramlist{i+1};
	ind = strmatch(lower(pname),optionsnames);
	if isempty(ind)
		error(['Invalid parameter: ''' pname '''.']);
	elseif length(ind) > 1
		error(['Ambiguous parameter: ''' pname '''.']);
	end
	switch(optionsnames{ind})
		case 'datafile'
			if ischar(pvalue)
				options.dataFile = pvalue;
			else
				error(sprintf(msgInvalidPair,pname));
            end
%             if ~isempty(findstr(pvalue , '.txt'))
% 			    options.dataType = 'text';
%             elseif ~isempty(findstr(pvalue, '.mat'))
%                 options.dataType = 'matlab';
%             elseif ~isempty(findstr(pvalue, '.xls'))
%                 options.dataType = 'excel';
%             else
%                 error('*** ERROR: unrecognized data format');
%             end
        case 'mixturetype'
            if ischar(pvalue)
                if ~strmatch(pvalue, strvcat('mix','ad_mix','linear_mix','codon_mix','spatical'),'exact')
                    error('*** ERROR: unrecoganized model type');
                end
                if isempty(pvalue),
                    options.mixtureType = '.';
                else
                    options.mixtureType = pvalue;
                end
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'initialk'
            pvalue = str2num(pvalue);
            if isnumeric(pvalue)
                if isempty(pvalue),
                    options.initialK = 0;
                else
                    options.initialK = pvalue;
                end
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'fixedk'
            if ischar(pvalue)
                if isempty(pvalue),
                    options.fixedK = 'no';
                else
                    options.fixedK = pvalue;
                end
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'namefile'
            if ischar(pvalue)
                options.nameFile = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'indexfile'
            if ischar(pvalue)
                options.indexFile = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'outputmat'
            if ischar(pvalue)
                options.outputMat = pvalue;
                directoryName = fileparts(pvalue);
                if ~exist(directoryName)
                   fprintf(1,'*** ERROR: Output directory ''%s'' does not exist.\n', directoryName);
                   return
                end
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'datatype'
            if ischar(pvalue)
                options.dataType = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'linkagemap'
            if ischar(pvalue)
                options.linkageMap = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'coordinatefile'
            if ischar(pvalue)
                options.coordinateFile = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'groups'
            if ischar(pvalue)
                options.groups = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'groupname'
            if ischar(pvalue)
                options.groupname = pvalue;
            else
                error(sprintf(msgInvalidPair,pname));
            end
            
        % the options below are for admixture analysis
        case 'clusters'
            if ischar(pvalue)
                options.clusters = str2num(pvalue);                    
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'minsize'
            if ischar(pvalue)
                options.minSize = str2num(pvalue);                    
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'iters'
            if ischar(pvalue)
                options.iters = str2num(pvalue);                    
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'refinds'
            if ischar(pvalue)
                options.refInds = str2num(pvalue);                    
            else
                error(sprintf(msgInvalidPair,pname));
            end
        case 'refiters'
            if ischar(pvalue)
                options.refIters = str2num(pvalue);                    
            else
                error(sprintf(msgInvalidPair,pname));
            end
		otherwise
			error(['Invalid parameter: ''' pname '''.']);
	end
end

% The subfunction to check syntax
if ~checkSyntax(options)
    return
end

switch options.mixtureType
    case 'mix'
        if isequal(options.groups,'yes')
            greedyPopMix_parallel(options);
        else
            independent_parallel(options);
        end
    case 'linear_mix'
        linkage_parallel(options);
    case 'codon_mix'
        linkage_parallel(options);
    case 'spatial'
        if isequal(options.groups, 'yes')
            spatialPopMixture_parallel(options);
        else
            spatial_parallel(options);
        end
    case 'ad_mix'
        admix_parallel(options);
end

% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------
function isOK = checkSyntax(options)
isOK = 1;
if strcmp(options.fixedK, 'yes') && length(options.initialK)>1
    display('*** ERROR: conflicting in options fixedk and initialk.');
    isOK = 0;
end

if strcmp(options.mixtureType, 'mix')
    if strcmp(options.dataType, 'excel') || strcmp(options.dataType,'sequence')
        display('*** ERROR: unknown datatype for the independence module.');
        isOK = 0;
    end
end

% check the admixture parameters
admix_str = {options.clusters, options.minSize, options.iters, ...
                                        options.refInds, options.refIters};
pt = cellfun('isempty', admix_str);

if all(pt)
    if ~strcmp(options.mixtureType, 'ad_mix')
        isOK = 1;
    else
        display('*** ERROR: problematic mixture type.');
        isOK = 0;
    end
end

if any(pt) && strcmp(options.mixtureType, 'ad_mix')
    display('*** ERROR: incomplete admixture parameters.');
    isOK = 0;
else
    isOK = 1;
end

