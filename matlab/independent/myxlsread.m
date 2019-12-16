function [A B] = myxlsread(file)
% This function read a tab ('\t') seperated txt file
% input file structure:
% first row: title
% sencond to end row: first column, sample ID
%                     second column, cluster label
%                     other columns, gene sequences
% Lu Cheng
% 26.06.2010

% there can be multiple numeric columns in the input file
% Lu Cheng, 25.11.2010

delimiter = '\t';

if exist(file,'file')~=2
    error('The input file %s does not exist!', file);
end

lines = textread(file,'%s','delimiter','\n');

title = strread(lines{1},'%s','delimiter',delimiter);
nRow = length(lines);
nCol = length(title);

% determine numeric Columns
tmp = strread(lines{2},'%s','delimiter',delimiter);
numCols = [];
for i = 1:length(tmp)
    if ~isnan(str2double(tmp{i}))
        numCols(end+1) = i;  %#ok<AGROW>
    end
end

A = cell(nRow-1, length(numCols));
B = cell(nRow, nCol);

B(1,:) = title;
for i=2:nRow
    if isempty(lines{i})
        B(i,:) = [];
        A(i-1,:) = [];
    else
        B(i,:) = strread(lines{i},'%s','delimiter',delimiter);
        A(i-1,:) = B(i,numCols);
    end
end

A = cellfun(@str2double,A);