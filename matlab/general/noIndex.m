function puredata = noIndex(data, noalle)
% NOINDEX Check that the data contains no index column.
%  Input: two variables from a mixture/admixture result structure.
%  Output: 
%        puredata: a data contains no index column.

if size(data,2) == length(noalle) + 1
    puredata = data(:,[1:end-1]); % remove the index column
else
    puredata = data;
end
