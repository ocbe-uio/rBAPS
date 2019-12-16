function iconn
% Change the icon of all open figures to match the figure number
%
% iconn
%
% See also: seticon
%
% In order to have this feature automatically on for all new figures
% execute the following command:
%   set(0,'defaultfigurecreatefcn','iconn')
% You may want to insert it into your startup.m file

% Copyright 2000-2002, Research and Development

h = get(0,'children');
for i=1:min(length(h), 9)
    seticon(i, which(sprintf('icon%d.ico', i)))
end
