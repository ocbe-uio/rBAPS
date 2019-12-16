function [NT,AA] = seqcode()
%SEQCODE - Return vector for mapping sequence letters to integers
%
% Syntax: [NT,AA] = seqcode
%
% See also:

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

NT = 'ACGTDI?-';
if (nargout>1)
AA = 'ARNDCQEGHILKMFPSTWYV*-';
end

% AANames = {'ala' 'arg' 'asn' 'asp' 'cys' 'gln' 'glu' 'gly' 'his' 'ile' 'leu' 'lys' 'met' ...
%            'phe' 'pro' 'ser' 'thr' 'trp' 'tyr' 'val'};