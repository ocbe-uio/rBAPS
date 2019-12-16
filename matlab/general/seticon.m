function seticon(window, useicon)
% Set icon for window.
%
% seticon(window, useicon)
%
% Parameters:
%  window: Figure number or name of window
%  useicon: Icon number to use or file name of icon
%           =1 : Application
%           =2 : Hand        (x)
%           =3 : Question    (?)
%           =4 : Exclamation /!\
%           =5 : Asterisk    (i)
%           =6 : Winlogo
%
% Examples:
%  seticon(2,3)
%  seticon(1, 'iconfile.ico')
%  seticon('Microsoft Internet', 6)

% Copyright 2000-2002, Research and Development

if nargin~=2
    error('Two arguments required')
end

if ~isstr(window)
    window = wgetname(window);
end

if ~any(window)
    warning('Window specification insufficient')
    return
end

if isstr(useicon)
    switch lower(useicon)
    case 'application'
        icon(1, window);
    case 'hand'
        icon(2, window);
    case 'question'
        icon(3, window);
    case 'exclamation'
        icon(4, window);
    case 'asterisk'
        icon(5, window);
    case 'winlogo'
        icon(6, window);
    otherwise
        icon(101, window, useicon);
    end
else
    if useicon>=1 & useicon<=6
        icon(useicon, window);
    else
        error('Icon number out of range')
    end
end
