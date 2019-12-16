function display(tr)
%DISPLAY command window display of phylogenetic tree objects.

% Copyright 2003-2005 The MathWorks, Inc.
% $Revision: 1.1.6.4 $ $Author: batserve $ $Date: 2005/06/09 21:55:53 $

n = numel(tr);
if n > 1
    disp(tr)
elseif n==0
    disp('    Empty array of phylogenetic tree objects')
else
    n = length(tr.dist); 
    switch n
        case 0
            disp('    Empty phylogenetic tree object')
        case 1
            disp('    Phylogenetic tree object with 1 leaf (0 branches)')
        case 3
            disp('    Phylogenetic tree object with 2 leaves (1 branch)')
        otherwise
            disp(['    Phylogenetic tree object with ' num2str((n+1)/2) ...
                  ' leaves (' num2str((n-1)/2) ' branches)'])
    end
end

