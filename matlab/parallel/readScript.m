function readScript(filename)
% READSCRIPT read the script file and output the parameters
% this function does not perform syntax checking.
% Example:
% readScript('script.txt')

% read the script
ind = readfile(filename);
if isempty(ind)
    return
end
nLines = size(ind,1);

% extract command information
optionStr = [];
for k = 1:nLines
    [cmdName, paraStr] = extract(ind(k,:));
    optionStr = [optionStr cmdName ',' paraStr ','];
end
optionStr = optionStr(1:end-1); % remove the last coma

% call function parallel
eval(['parallel(' optionStr ')'])

% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------
function [cmdName, paraStr] = extract(commandline)
% function to extract the command name and the parameter string

[cmdName, remainStr] = strtok(commandline,'(');
boundary = regexp(remainStr,'''');

if isempty(boundary) % if paraStr does not contain quotation marks
    % use parenthesis as boundaries
    startPt = regexp(remainStr,'(') + 1;
    endPt = regexp(remainStr,')') - 1;
else
    startPt = boundary(1) + 1;
    endPt = boundary(2) - 1;
end
paraStr = remainStr(startPt: endPt);

cmdName = strcat('''',cmdName,'''');
paraStr = strcat('''',paraStr,'''');

% -------------------------------------------------------------------------
function T = readfile(filename);
f = fopen(filename,'r');
if f == -1 
    % error(filename); 
    display('*** ERROR: invalid script name.');
    T = [];
    return
end
i = 1;
while 1
  clear line;
  line = fgetl(f);
  if ~isstr(line), break, end  
  n = length(line);
  T(i,1:n) = line(1:n);
  i = i+1;
end
fclose(f);

