function S=i_encode_n(Seq)
%I_ENCODE_N - Convert a nucleotide sequence from a letter to an integer representation
%Internal function encodes the nucleotide sequences by digits
%
% Syntax: S=i_encode_n(Seq)
%
% Inputs:
%    Seq   - Letter representation of sequence
%
% Outputs:
%    S     - Integer representation of sequence
%
%
% See also: I_ENCODE_A

% Molecular Biology & Evolution Toolbox, (C) 2005
% Author: James J. Cai
% Email: jamescai@hkusua.hku.hk
% Website: http://web.hku.hk/~jamescai/
% Last revision: 5/28/2005
% seqcode = 'ACGTDI?-'
method=1;
Seq=upper(Seq); % Read lower case letters
[NT,AA] = seqcode; % NT = 'ACGTDI?-'
switch (method)
    case (1)
        [n,m]=size(Seq);
        S = ones(n,m).*8;
        Seq(find(Seq=='U'))='T';    %replace U with T
        for (k=1:8),
            S(find(Seq==NT(k)))=k;
        end
        S(find(Seq==NT(7))) = -999; % missing data denoted as -999.
        S(find(Seq==NT(8))) = -999; % gap 
    case (2)
        [n,m]=size(Seq);
        S = zeros(n,m);  i = 1; j = 1;
        
        code = zeros(256,1); 
        for o = 1:256
            code(o) = nan;
        end
        
        % NT = 'ACGT-';
        
        for o = 1:5
            code(abs(NT(o))) = o;  
        end
        
        for i=1:n
            for j=1:m
                if Seq(i,j) == 'U' 
                    S(i,j) = code(abs('T'));
                else
                    S(i,j) = code(abs(Seq(i,j)));
                end;
            end;
        end;
end

% S=uint8(S);