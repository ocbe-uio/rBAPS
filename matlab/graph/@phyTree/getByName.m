function sel = getbyname(tr,query,varargin)
%GETBYNAME Selects branches and leaves by name.
%
%   S = GETBYNAME(T,EXPRESSION) returns a logical vector S of size 
%   [NUMNODES x 1] indicating the node names of the phylogenetic tree T
%   that match the regular expression EXPRESSION regardless of case. 
%
%   Symbols than can be used in a matching regular expression are explained
%   in help REGEXP.
%
%   When EXPRESSION is a cell array of strings, GETBYNAME returns a matrix
%   where every column corresponds to every query in EXPRESSION.
%
%   S = GETBYNAME(T,STRING,'EXACT',true) looks for exact matches only
%   (ignoring case). When STRING is a cell array of strings, GETBYNAME
%   returns a vector with indices.
%
%   Example:
%
%      % Load a phylogenetic tree created from a protein family:
%      tr = phytreeread('pf00002.tree');
%       
%      % Select all the 'mouse' and 'human' proteins:
%      sel = getbyname(tr,{'mouse','human'});
%      view(tr,any(sel,2));
%       
%   See also PHYTREE, PHYTREE/PRUNE, PHYTREE/SELECT, PHYTREE/GET.

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.6.5 $ $Author: batserve $ $Date: 2005/06/09 21:55:55 $

if numel(tr)~=1
     error('Bioinfo:phytree:getbyname:NoMultielementArrays',...
           'Phylogenetic tree must be an 1-by-1 object.');
end

doExactMatch = false;

if  nargin > 2
    okargs = {'exact',''};
    for j=1:2:nargin-2
        pname = varargin{j};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Bioinfo:phytree:getbyname:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:phytree:getbyname:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1
                   if nargin == 3
                       doExactMatch = true;
                   else
                       doExactMatch = opttf(varargin{j+1});
                       if isempty(doExactMatch)
                           error('Bioinfo:phytree:getbyname:InputOptionNotLogical',...
                                '%s must be a logical value, true or false.',...
                                 upper(char(okargs(k))));
                       end
                   end
            end %switch
        end %if
    end %for
end %if


numLabels = numel(tr.names);
if iscell(query)
    if doExactMatch
        sel = zeros(numLabels,1);
    else
        sel = false(numLabels,numel(query));
    end
    for ind = 1:numel(query)
        if doExactMatch
            sel(strcmpi(query{ind},tr.names)) = ind;
        else
            try
                regexpiOutput = regexpi(tr(:).names,query{ind});
            catch
                error('Bioinfo:phytree:getbyname:IncorrectRegularExpression',...
                    ['The query expression produced the following error in ' ...
                    'REGEXPI: \n%s'],lasterr);
            end
            sel(:,ind) = ~cellfun('isempty',regexpiOutput);
        end
    end
else % must be a single string of chars
    if doExactMatch
        sel = strcmpi(query,tr.names);
    else
        try
            regexpiOutput = regexpi(tr(:).names,query);
        catch
            error('Bioinfo:phytree:getbyname:IncorrectRegularExpression',...
                ['The query expression produced the following error in ' ...
                'REGEXPI: \n%s'],lasterr);
        end
        sel = ~cellfun('isempty',regexpiOutput);
    end
end
