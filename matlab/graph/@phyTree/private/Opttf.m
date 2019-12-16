function tf = opttf(pval)
%OPTTF determines whether input options are true or false

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.3.4.2 $   $Date: 2004/12/24 20:42:39 $


if islogical(pval)
    tf = all(pval);
    return
end
if isnumeric(pval)
    tf = all(pval~=0);
    return
end
if ischar(pval)
    truevals = {'true','yes','on','t'};
    k = any(strcmpi(pval,truevals));
    if k
        tf = true;
        return
    end
    falsevals = {'false','no','off','f'};
    k = any(strcmpi(pval,falsevals));
    if k
        tf = false;
        return
    end
end
% return empty if unknown value
tf = logical([]);