function openvar(name, tr) %#ok
%OPENVAR Opens a phylogenetic tree object for graphical editing.

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.1.6.3 $ $Author: batserve $ $Date: 2006/06/16 20:06:44 $

try
    view(tr);
catch
    % rethrows the error into a dlg window
    errordlg(lasterr, 'Inspection error', 'modal');
end
