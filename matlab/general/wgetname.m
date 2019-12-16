function s=wgetname(h)
% Returns the window name of a window
%
% s=wgetname(h)
%
% Parameters:
%  h : window number
%
% See also: seticon

% Copyright 2000-2002, Research and Development

if strcmp(get(h,'numbertitle'), 'on')
  if get(h,'name')
    s = sprintf('Figure %d: %s', h, get(h,'name'));
  else
    s = sprintf('Figure %d', h);
  end
else
  s = get(h,'name');
end
