function [D]=chooseDistance(aln)
% CHOOSEDISTANCE Display a quest dialog for choosing a genetic distance
%  input: a structure containing part of variables from the mixture result.
%  ouput: a symmetric distance matrix

if (isempty(aln))
    error('Need input.');
end

if aln.npops == 1
    disp('*** ERROR: no population structure is found.');
    return
end

ButtonName=questdlg('Choose a type of genetic distance?', ...
    'Select genetic distance', ...
    'KL','Nei','Hamming','KL');

switch ButtonName,
    case 'KL'
        disp('Using KL distance.');
        D = dn_kl(aln);
    case 'Nei'
        disp('Using Nei distance.');
        D = dn_nei(aln);
    case 'Hamming'
        disp('Using Hamming distance. Please wait...');
        D = dn_hamming(aln);
%     case 'LogDet'
%         disp('Using LogDet distance');
%         D= dn_logdet(aln);
    otherwise
        D = [];
end

D = D + D'; % make it symmetric
%--------------------------------------------------------------------------
% SUBFUNCTIONS
%--------------------------------------------------------------------------
function [dist_mat] = dn_kl(aln)
npops = aln.npops;
COUNTS = aln.COUNTS;
adjprior = aln.adjprior;
% data = noIndex(aln.data, aln.noalle);
% partition = aln.partition;

dist_mat = zeros(npops, npops);
maxnoalle = size(COUNTS,1);
nloci = size(COUNTS,2);
d = zeros(maxnoalle, nloci, npops);

prior = adjprior;
prior(find(prior==1))=0;
nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yhtä alleelia.
prior(1,nollia)=1;
for pop1 = 1:npops
    d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
    %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
end

for pop1 = 1:npops
    % rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
    for pop2 = 1:pop1-1
        dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
        div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
        div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
        div = (div12+div21)/2;
        dist_mat(pop1,pop2) = div;
    end
end

%--------------------------------------------------------------------------
function [dist_mat] = dn_nei(aln)

npops = aln.npops;
COUNTS = aln.COUNTS;
% adjprior = aln.adjprior;
% data = aln.data;
% partition = aln.partition;

dist_mat = zeros(npops, npops);
maxnoalle = size(COUNTS,1);
nloci = size(COUNTS,2);
d = zeros(maxnoalle, nloci, npops);

for pop1 = 1:npops
    d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))),maxnoalle,1);
    %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
end
for pop1 = 1:npops
    % rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
    for pop2 = 1:pop1-1
        dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
        div1 = sum(sum(dist1.*dist2));
        div2 = sqrt(sum(sum(dist1.^2)))*sqrt(sum(sum(dist2.^2)));
        div = -log(div1/div2);
        dist_mat(pop1,pop2) = div;
    end
end

%--------------------------------------------------------------------------
function [dist_mat] = dn_hamming(aln)

npops = aln.npops;
data = noIndex(aln.data, aln.noalle);
partition = aln.partition;
dist_mat = zeros(npops, npops);
for pop1 = 1:npops
    for pop2 = 1:pop1-1
        dist_mat(pop1,pop2) = hamming_dist(data(logical(partition==pop1),[1:end-1]),...
            data(logical(partition==pop2),[1:end-1]));
    end
end

function dist = hamming_dist(data1,data2)
[length1,nloci] = size(data1);
length2 = size(data2,1);
dist1 = 0;
for i = 1:length1
    dist2 = 0;
    for j = 1:length2
        dist2 = dist2 + sum(data1(i,:)~=data2(j,:))/nloci;
    end
    dist1 = dist1 + dist2/length2;
end
dist = dist1/length1;

%--------------------------------------------------------------------------
function [D]=dn_logdet(aln)
%DN_LOGDET - Log-det (paralinear) distance
%The LogDet model computes the distance from the determinant of the matrix of
%co-occurrence of nucleotides in the two species, according to the formula
%
%   D  = - 1/4(loge(|F|) - 1/2loge(fA1 fC1 fG1 fT1 fA2 fC2 fG2 fT2))
%
%Where F is a matrix whose (i,j) element is the fraction of sites at which base
%i occurs in one species and base j occurs in the other. fji is the fraction of
%sites at which species i has base j. The LogDet distance cannot cope with
%ambiguity codes. It must have completely defined sequences. One limitation of
%the LogDet distance is that it may be infinite sometimes, if there are too many
%changes between certain pairs of nucleotides. This can be particularly
%noticeable with distances computed from bootstrapped sequences.
%
% Syntax: [D]=dn_logdet(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    D     - Distance matrix
%
% See also:

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 5/3/2006

if (isstruct(aln)),
    S=aln.data;
else
    S=aln;
end

[n,m] = size(S);
D = zeros(n,n);

for i=1:n-1
    for j=i+1:n
        D(i,j) = d_logdet(S(i,:), S(j,:));
        D(j,i) = D(i,j);
    end
end

%--------------------------------------------------------------------------
function d=d_logdet(seq1, seq2)
[S,gap] = countntchange(seq1, seq2);
% S=count_ntchanges(seq1, seq2);
if (det(S)<=0)
    d=inf;
else
    f1=sum(S);
    f2=sum(S');
    d=(-1/4)*(log(det(S))-(1/2)*log( prod(f1)*prod(f2)));
end

%--------------------------------------------------------------------------
function [D,gap]=countntchange(s1,s2)
%COUNTNTCHANGE - Count nucleotide changes in two DNA sequences
%D is a 4x4 array, with bases in seq1 along top, seq2 along side,
%in order A,C,G,T.
%
% Syntax: [D,gap]=countntchange(s1,s2)
%
% Inputs:
%    s1      - Sequence 1 vector
%    s2      - Sequence 2 vector
%
% Outputs:
%    D       - Codon Adaptation Index value
%    gap     - Codon Adaptation Index value
%
%
% See also: COUNTAACHANGE COUNTCDCHANGE

% Molecular Biology & Evolution Toolbox, (C) 2006
% Author: James J. Cai
% Email: jamescai@hku.hk
% Website: http://bioinformatics.org/mbetoolbox/
% Last revision: 3/28/2006

if (nargout>1)
    [D,gap]=countchange(s1,s2,4);
else
    [D]=countchange(s1,s2,4);
end

% -------------------------------------------------------------------------
function [D,gap] = countchange(s1,s2,nword)

if ~(ismember(nword, [4 20 61])),
    error('Wrong NWORD')
end
if (length(s1)~=length(s2)),
    error('Sequences are not of same length.')
end

D=zeros(nword);
WORD=1:nword;

for j=1:nword
    s1sites=(s1==WORD(j));
    for i=1:nword
        D(i,j)=sum(s1sites & (s2==WORD(i)));
    end
end

if (nargout>1),
    gap = length(s1)-sum(D(:));	% gaps
end

