function [aln2] = encodealn(aln)
%ENCODEALN - Convert nucleotide in alignment to integer.
%
% Syntax: [aln2] = encodealn(aln)
%
% Inputs:
%    aln     - Alignment structure letter representation
%
% Outputs:
%    aln2    - Alignment structure integer representation
%
%
% See also: CODONISESEQ, ENCODESEQ

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005

if ~(aln.seqtype) error('Do not know the type of sequence! ... BAPS'); end
aln2=aln;

switch (aln2.seqtype)
    case (1)
         aln2.seq = i_encode_n(aln2.seq);
    case (2)
         aln2.seq = i_encode_n(aln2.seq);
    case (3)
         aln2.seq = i_encode_a(aln2.seq);
end