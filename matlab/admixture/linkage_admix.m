function linkage_admix(tietue)

global COUNTS; global PARTITION; global SUMCOUNTS;
clearGlobalVars;

PARTITION = tietue.PARTITION;
COUNTS = tietue.COUNTS;
SUMCOUNTS = tietue.SUMCOUNTS;
rowsFromInd = tietue.rowsFromInd;
data = double(tietue.data);
npops = tietue.npops;
noalle = tietue.noalle;
switch tietue.mixtureType
    case 'linear_mix'
        linkage_model = 'linear';
    case 'codon_mix'
        linkage_model = 'codon';
end
if isfield(tietue,'gene_lengths')
    gene_lengths = tietue.gene_lengths;
else
    [filename, pathname] = uigetfile('*.txt', 'Load file with lengths of the genes (same order as in data).');
    gene_lengths = load([pathname filename]);
end

if length(unique(PARTITION(find(PARTITION>0))))==1
    disp('Only one population in the input file');
    disp('No admixture detected');
    return
end

answers = inputdlg({['Input the minimum size of a population that will'...
            ' be taken into account when admixture is estimated.']},...
            'Input minimum population size',1,{'5'});
if isempty(answers), return; end
alaRaja = str2double(answers{1,1});
[npops] = poistaLiianPienet(npops, rowsFromInd, alaRaja);

nloci = size(COUNTS,2);
ninds = size(data,1)/rowsFromInd;

answers = inputdlg({'Input number of iterations'},'Input',1,{'50'});
if isempty(answers), return; end
iterationCount = str2double(answers{1,1});

answers = inputdlg({'Input number of reference individuals from each population'},'Input',1,{'50'});
if isempty(answers), nrefIndsInPop = 50;
else nrefIndsInPop = str2double(answers{1,1});
end

answers = inputdlg({'Input number of iterations for reference individuals'},'Input',1,{'10'});
if isempty(answers), return; end
iterationCountRef = str2double(answers{1,1});

[cliq_data, sep_data, cliq_counts, component_mat] = createCliqData(data, gene_lengths, noalle, ...
    linkage_model, rowsFromInd);

% Repeat: simulate clique frequencies and estimate proportions. Save the
% average proportions in "proportionsIt".

proportionsIt = zeros(ninds,npops);
for iterationNum = 1:iterationCount
    disp(['Iter: ' num2str(iterationNum)]);
    %allfreqs = simulateAllFreqs(noalle);
    [cliq_freqs, sep_freqs] = simulateCliqFreqs(cliq_counts, noalle, component_mat, gene_lengths, linkage_model);
    for ind=1:ninds

        %omaFreqs = computePersonalAllFreqs(ind, data, allfreqs, rowsFromInd);
        [ownCliqFreqs, ownSepFreqs] = computePersonalCliqueFreqs(ind, ...
            cliq_data, cliq_freqs, sep_data, sep_freqs, rowsFromInd, gene_lengths, linkage_model);
        osuusTaulu = zeros(1,npops);
        if PARTITION(ind)==0
            % Outlier individual
        elseif PARTITION(ind)~=0
            if PARTITION(ind)>0
                osuusTaulu(PARTITION(ind)) = 1;
            else
                % Individuals who are not assigned to any cluster.
                arvot = zeros(1,npops);
                for q=1:npops
                    osuusTaulu = zeros(1,npops);
                    osuusTaulu(q) = 1;
                    arvot(q) = computeIndLikelihood(ownCliqFreqs, ownSepFreqs, osuusTaulu);
                end
                [iso_arvo, isoimman_indeksi] = max(arvot);
                osuusTaulu = zeros(1,npops);
                osuusTaulu(isoimman_indeksi) = 1;
                PARTITION(ind)=isoimman_indeksi;
            end
            logml = computeIndLikelihood(ownCliqFreqs, ownSepFreqs, osuusTaulu);
            
            for osuus = [0.5 0.25 0.05 0.01]
                [osuusTaulu, logml] = searchBest(osuus, osuusTaulu, ownCliqFreqs, ownSepFreqs, logml);
            end
        end
        proportionsIt(ind,:) = proportionsIt(ind,:).*(iterationNum-1) + osuusTaulu;
        proportionsIt(ind,:) = proportionsIt(ind,:)./iterationNum;
    end
end

disp(['Creating ' num2str(nrefIndsInPop) ' reference individuals from ']);
disp('each population.');

%allfreqs = simulateAllFreqs(noalle);  % Simuloidaan alleelifrekvenssisetti
%allfreqs = computeAllFreqs2(noalle); % Koitetaan tällaista.
%refData = simulateIndividuals(nrefIndsInPop,rowsFromInd,allfreqs);

exp_cliq_freqs = ...
    computeExpectedFreqs(cliq_counts, noalle, component_mat, gene_lengths, linkage_model);
[ref_cliq_data, ref_sep_data] = ...
    simulateLinkageIndividuals(nrefIndsInPop, rowsFromInd, exp_cliq_freqs, ...
    gene_lengths, noalle, component_mat, linkage_model);
nrefInds = npops*nrefIndsInPop;

disp(['Analysing the reference individuals in ' num2str(iterationCountRef) ' iterations.']);
refProportions = zeros(nrefInds,npops);
for iter = 1:iterationCountRef
    disp(['Iter: ' num2str(iter)]);
    %allfreqs = simulateAllFreqs(noalle);
    [cliq_freqs, sep_freqs] = simulateCliqFreqs(cliq_counts, noalle, component_mat, gene_lengths, linkage_model);
    for ind = 1:nrefInds
        %omaFreqs = computePersonalAllFreqs(ind, refData, allfreqs, rowsFromInd);
        [ownCliqFreqs, ownSepFreqs] = computePersonalCliqueFreqs(ind, ...
            ref_cliq_data, cliq_freqs, ref_sep_data, sep_freqs, rowsFromInd, gene_lengths, linkage_model);
        osuusTaulu = zeros(1,npops);
        pop = ceil(ind/nrefIndsInPop);
        osuusTaulu(pop)=1;
        %logml = computeIndLogml(omaFreqs, osuusTaulu);
        logml = computeIndLikelihood(ownCliqFreqs, ownSepFreqs, osuusTaulu);
        for osuus = [0.5 0.25 0.05 0.01]
            %[osuusTaulu, logml] = searchBest(osuus, osuusTaulu, omaFreqs, logml);
            [osuusTaulu, logml] = searchBest(osuus, osuusTaulu, ownCliqFreqs, ownSepFreqs, logml);
        end
        refProportions(ind,:) = refProportions(ind,:).*(iter-1) + osuusTaulu;
        refProportions(ind,:) = refProportions(ind,:)./iter;
    end
end
refTaulu = zeros(npops,100);
for ind = 1:nrefInds
    pop = ceil(ind/nrefIndsInPop);
    omanOsuus = refProportions(ind,pop);
    if round(omanOsuus*100)==0
        omanOsuus = 0.01;
    end
    if abs(omanOsuus)<1e-5
        omanOsuus = 0.01;
    end
    refTaulu(pop, round(omanOsuus*100)) = refTaulu(pop, round(omanOsuus*100))+1;
end
    
% Rounding:
proportionsIt = proportionsIt.*100; proportionsIt = round(proportionsIt);
proportionsIt = proportionsIt./100;
for ind = 1:ninds
    % if sum not equal to one, fix the largest part.
    if (PARTITION(ind)>0) && (sum(proportionsIt(ind,:)) ~= 1)
        [isoin,indeksi] = max(proportionsIt(ind,:));
        erotus = sum(proportionsIt(ind,:))-1;
        proportionsIt(ind,indeksi) = isoin-erotus;
    end
end

% "p-value" for admixture
uskottavuus = zeros(ninds,1);
for ind = 1:ninds
    pop = PARTITION(ind);
    if pop==0  % an outlier
        uskottavuus(ind)=1;
    else
        omanOsuus = proportionsIt(ind,pop);
        if abs(omanOsuus)<1e-5
            omanOsuus = 0.01;
        end
        if round(omanOsuus*100)==0
            omanOsuus = 0.01;
        end
        refPienempia = sum(refTaulu(pop, 1:round(100*omanOsuus)));
        uskottavuus(ind) = refPienempia / nrefIndsInPop;
    end
end

tulostaAdmixtureTiedot(proportionsIt, uskottavuus, alaRaja, iterationCount); 

%viewPartition(proportionsIt, popnames);

[filename, pathname] = uiputfile('*.mat','Save admixture results as');
if (filename == 0) & (pathname == 0)
    % Cancel was pressed
    return
end

% copy 'baps4_output.baps' into the text file with the same name.
if exist('baps4_output.baps','file')
    copyfile('baps4_output.baps',[pathname filename '.txt'])
    delete('baps4_output.baps')
end

tietue.proportionsIt = proportionsIt;
tietue.pvalue = uskottavuus; % Added by Jing
tietue.admixnpops = npops;
tietue.mixtureType = 'admix'; % added by jing on 09.09.2008
% save([pathname filename], 'tietue');
save([pathname filename], 'tietue','-v7.3'); % added by Lu Cheng, 08.06.2012


%----------------------------------------------------------------------------


function [npops] = poistaLiianPienet(npops, rowsFromInd, alaraja)
% Muokkaa tulokset muotoon, jossa outlier yksilät on
% poistettu. Tarkalleen ottaen poistaa ne populaatiot, 
% joissa on vähemmän kuin 'alaraja':n verran yksiläit?

global PARTITION;
global COUNTS;
global SUMCOUNTS;

popSize=zeros(1,npops);
for i=1:npops
    popSize(i)=length(find(PARTITION==i));
end
miniPops = find(popSize<alaraja);

if length(miniPops)==0
    return;
end

outliers = [];
for pop = miniPops
    inds = find(PARTITION==pop);
    disp('Removed individuals: ');
    disp(num2str(inds));
    outliers = [outliers; inds];
end

ninds = length(PARTITION);
PARTITION(outliers) = 0;
korit = unique(PARTITION(find(PARTITION>0)));
for n=1:length(korit)
    kori = korit(n);
    yksilot = find(PARTITION==kori);
    PARTITION(yksilot) = n;
end
COUNTS(:,:,miniPops) = [];
SUMCOUNTS(miniPops,:) = [];

npops = npops-length(miniPops);

%------------------------------------------------------------------------

function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];

%--------------------------------------------------------


function allFreqs = computeAllFreqs2(noalle)
% Lisää a priori jokaista alleelia
% joka populaation joka lokukseen j 1/noalle(j) verran.

global COUNTS;
global SUMCOUNTS;

max_noalle = size(COUNTS,1);
nloci = size(COUNTS,2);
npops = size(COUNTS,3);

sumCounts = SUMCOUNTS+ones(size(SUMCOUNTS));
sumCounts = reshape(sumCounts', [1, nloci, npops]);
sumCounts = repmat(sumCounts, [max_noalle, 1 1]);

prioriAlleelit = zeros(max_noalle,nloci);
for j=1:nloci
    prioriAlleelit(1:noalle(j),j) = 1/noalle(j);
end
prioriAlleelit = repmat(prioriAlleelit, [1,1,npops]);
counts = COUNTS + prioriAlleelit;
allFreqs = counts./sumCounts;

%--------------------------------------------------------------------------


function refData = simulateIndividuals(n,rowsFromInd,allfreqs)
% simuloidaan n yksilää jokaisesta populaatiosta. (

npops = size(allfreqs,3);
nloci = size(allfreqs,2);
ninds = n*npops;

refData = zeros(ninds*rowsFromInd,nloci);
counter = 1;  % Pitää kirjaa mille riville seuraavaksi simuloidaan.

for ind = 1:ninds
    pop = ceil(ind/n);
    for loc = 1:nloci
        for k=0:rowsFromInd-1
            refData(counter+k,loc) = simuloiAlleeli(allfreqs,pop,loc);
        end
    end
    counter = counter+rowsFromInd;
end

function all = simuloiAlleeli(allfreqs,pop,loc)
% Simuloi populaation pop lokukseen loc alleelin.
freqs = allfreqs(:,loc,pop);
cumsumma = cumsum(freqs);
arvo = rand;
isommat = find(cumsumma>arvo);
all = min(isommat);

%---------------------------------------------------------------------------


function loggis = computeIndLikelihood(ownCliqFreqs, ownSepFreqs, proportions)
% Calculates the likelihood when the origins are defined by "proportions".

aux = proportions * ownCliqFreqs;
aux = log(aux);
loggis = sum(aux);

clear aux;

aux2 = proportions * ownSepFreqs;
aux2 = log(aux2);
loggis = loggis - sum(aux2);


%--------------------------------------------------------------------------


function osuusTaulu = suoritaMuutos(osuusTaulu, osuus, indeksi)
% Päivittää osuusTaulun muutoksen jälkeen.

global COUNTS;
npops = size(COUNTS,3);

i1 = rem(indeksi,npops);
if i1==0, i1 = npops; end;
i2 = ceil(indeksi / npops);

osuusTaulu(i1) = osuusTaulu(i1)-osuus;
osuusTaulu(i2) = osuusTaulu(i2)+osuus;


%-------------------------------------------------------------------------


function [osuusTaulu, logml] = searchBest(osuus, osuusTaulu, ownCliqFreqs, ownSepFreqs, logml)

ready = 0;
while ready ~= 1
    muutokset = calcChanges(osuus, osuusTaulu, ownCliqFreqs, ownSepFreqs, logml);
    [maxMuutos, indeksi] = max(muutokset(1:end));
    if maxMuutos>0
        osuusTaulu = suoritaMuutos(osuusTaulu, osuus, indeksi);
        logml = logml + maxMuutos;
    else
        ready = 1;
    end
end



%---------------------------------------------------------------------------


function muutokset = calcChanges(osuus, osuusTaulu, ownCliqFreqs, ownSepFreqs, logml)
% Palauttaa npops*npops taulun, jonka alkio (i,j) kertoo, mik?on
% muutos logml:ss? mikäli populaatiosta i siirretään osuuden verran
% todennäkäisyysmassaa populaatioon j. Mikäli populaatiossa i ei ole
% mitään siirrettävää, on vastaavassa kohdassa rivi nollia.

global COUNTS;
npops = size(COUNTS,3);

notEmpty = find(osuusTaulu>0.005);
muutokset = zeros(npops);
empties = ~notEmpty;

for i1=notEmpty
    
    osuusTaulu(i1) = osuusTaulu(i1)-osuus;
    
    for i2 = [1:i1-1 i1+1:npops]
        osuusTaulu(i2) = osuusTaulu(i2)+osuus;
        loggis = computeIndLikelihood(ownCliqFreqs, ownSepFreqs, osuusTaulu);
        muutokset(i1,i2) = loggis-logml;
        osuusTaulu(i2) = osuusTaulu(i2)-osuus; 
    end
    
    osuusTaulu(i1) = osuusTaulu(i1)+osuus;
end


%---------------------------------------------------------------


function dispLine
disp('---------------------------------------------------');


%--------------------------------------------------------------------------


function tulostaAdmixtureTiedot(proportions, uskottavuus, alaRaja, niter)
h0 = findobj('Tag','filename1_text');
inputf = get(h0,'String');
h0 = findobj('Tag','filename2_text');
outf = get(h0,'String'); clear h0;

if length(outf)>0
    fid = fopen(outf,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
end

ninds = length(uskottavuus);
npops = size(proportions,2);
disp(' ');
dispLine;
disp('RESULTS OF ADMIXTURE ANALYSIS BASED');
disp('ON MIXTURE CLUSTERING OF INDIVIDUALS');
disp(['Data file: ' inputf]);
disp(['Number of individuals: ' num2str(ninds)]);
disp(['Results based on ' num2str(niter) ' simulations from posterior allele frequencies.']);
disp(' ');
if fid ~= -1
    fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['--------------------------------------------']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['RESULTS OF ADMIXTURE ANALYSIS BASED']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['ON MIXTURE CLUSTERING OF INDIVIDUALS']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Data file: ' inputf]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Number of individuals: ' num2str(ninds)]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Results based on ' num2str(niter) ' simulations from posterior allele frequencies.']); fprintf(fid, '\n');
    fprintf(fid, '\n');
end

ekaRivi = blanks(6);
for pop = 1:npops
    ekaRivi = [ekaRivi blanks(3-floor(log10(pop))) num2str(pop) blanks(2)];
end
ekaRivi = [ekaRivi blanks(1) 'p']; % Added on 29.08.06
disp(ekaRivi);
for ind = 1:ninds
    rivi = [num2str(ind) ':' blanks(4-floor(log10(ind)))];
    if any(proportions(ind,:)>0)
        for pop = 1:npops-1
            rivi = [rivi proportion2str(proportions(ind,pop)) blanks(2)];
        end
        rivi = [rivi proportion2str(proportions(ind,npops)) ':  '];
        rivi = [rivi ownNum2Str(uskottavuus(ind))];
    end
    disp(rivi);
    if fid ~= -1
        fprintf(fid,'%s \n',[rivi]); fprintf(fid,'\n');
    end
end
if fid ~= -1
    fclose(fid);
else
    diary off
end

%------------------------------------------------------

function str = proportion2str(prob)
%prob belongs to [0.00, 0.01, ... ,1]. 
%str is a 4-mark presentation of proportion.

if abs(prob)<1e-3
    str = '0.00';
elseif abs(prob-1) < 1e-3;
    str = '1.00';
else
    prob = round(100*prob);
    if prob<10
        str = ['0.0' num2str(prob)];    
    else
        str = ['0.' num2str(prob)];
    end;        
end;

%-------------------------------------------------------

function g=randga(a,b)
flag = 0;
if a>1 
c1 = a-1; c2 = (a-(1/(6*a)))/c1; c3 = 2/c1; c4 = c3+2; c5 = 1/sqrt(a); 
U1=-1; 
while flag == 0, 
if a<=2.5, 
U1=rand;U2=rand; 
else 
while ~(U1>0 & U1<1), 
U1=rand;U2=rand; 
U1 = U2 + c5*(1-1.86*U1); 
end %while 
end %if
W = c2*U2/U1; 
if c3*U1+W+(1/W)<=c4, 
flag = 1; 
g = c1*W/b; 
elseif c3*log(U1)-log(W)+W<1, 
flag = 1; 
g = c1*W/b; 
else 
U1=-1; 
end %if 
end %while flag
elseif a==1 
g=sum(-(1/b)*log(rand(a,1))); 
else 
while flag == 0, 
U = rand(2,1);
if U(1)>exp(1)/(a+exp(1)), 
g = -log(((a+exp(1))*(1-U(1)))/(a*exp(1))); 
if U(2)<=g^(a-1), 
flag = 1; 
end %if 
else 
g = ((a+exp(1))*U(1)/((exp(1))^(1/a))); 
if U(2)<=exp(-g), 
flag = 1; 
end %if 
end %if 
end %while flag 
g=g/b; 
end %if;


%-------------------------------------------------

function svar=randdir(counts,nc)
% Käyttäesim randdir([10;30;60],3)

svar=zeros(nc,1);
for i=1:nc
   svar(i,1)=randga(counts(i,1),1);
end
svar=svar/sum(svar);

%--------------------------------------------------

function waitALittle
A = rand(500);
gammaln(A);


%--------------------------------------------------

function [cliq_data, sep_data, cliq_counts, component_mat] = ...
    createCliqData(data, gene_lengths, noalle, linkage_model, ...
    rowsFromInd)
% cliq_data: cell array, each cell corresponds to one gene. Element (i,j)
% in cell k is the code of the allele combination in the j:th clique in
% gene k, for individual i.

% sep_data: like cliq_data. i:th separator separates cliques i and i+1.

% cliq_counts: cell array, each cell corresponds to one gene. Each cell is
% a 3-dimensional array, where the element (i,j,k) is the observed count of
% allele combination i, in clique j, in population k.

% The coding of the allele combinations: If a clique of 3 sites has
% noalle values 3,2,4, then the allele combinations are given numbers in
% lexicographic order: 111, 112, 113, 114, 121, 122, ..., 324.

%------------------------------------------------------------------------

% cliq_data on cell-array, jossa kukin solu vastaa yht?geeni? Alkio (i,j)
% solussa k merkitsee sen alleelikombinaation koodia, joka yksiläll?i
% havaitaan geenin k klikiss?numero j.

% cliq_counts on cell-array, jossa myäs kukin solu vastaa yht?geeni?
% Kukin solu on kolmiulotteinen taulukko, jonka alkio (i,j,k) on
% populaatiossa k, ko geenin j:nness?klikiss?havaitun alleelikombinaation
% i lukumäär?

% Alleelikombinaatioiden koodaus: Jos kolmen position klikiss?on (koko
% datassa) noalle:t 3,2,4, (eli ekassa positiossa alleelit 1-3, tokassa 1-2
% ja kolmannessa alleelit 1-4), niin alleelikombinaatiot numeroidaan
% leksikograafisessa järjestyksess? 111, 112, 113, 114, 121, 122, ...,
% 324.

%-----------------------------------------------------------------------

global PARTITION;

if sum(gene_lengths) ~= size(data,2)
    disp('Error 155');
end
if ~isa(data,'double')
    data = double(data);  % Required in matlab 6
end

ninds = size(data,1);
n_genes = length(gene_lengths);
if strcmp(linkage_model,'codon')
    n_cliques = gene_lengths-2;  % Number of cliques in each gene
else
    n_cliques = gene_lengths-1;  % Use linear model
end
max_noalle = zeros(n_genes,1);  % Maximum "clique noalle" in each gene.

component_mat = zeros(n_genes, max(gene_lengths));
cum_length = cumsum(gene_lengths);
component_mat(1,1:gene_lengths(1))=1:gene_lengths(1);
for i = 2:n_genes
    component_mat(i,1:gene_lengths(i)) = cum_length(i-1)+1:cum_length(i);
end

for i = 1:n_genes
    % What is the largest number of different values that are observed for
    % some clique in this gene:
    number = 0;
    if n_cliques(i)<1
        % the gene is shorter than a normal clique.
        if gene_lengths(i)==2
            % linkage_model must be 'codon' to end up here..
            positions = component_mat(i, [1 2]);
            number = prod(noalle(positions));
        else
            % gene_lengths(i) == 1
            positions = component_mat(i,1);
            number = noalle(positions);
        end
    else
        for j = 1:n_cliques(i)
            if strcmp(linkage_model,'codon'), positions = component_mat(i , j:j+2);
            else positions = component_mat(i, j:j+1);
            end
            
            cand = prod(noalle(positions));
            if cand>number
                number=cand;
            end
        end
    end
    max_noalle(i) = number;
end
        
cliq_data = cell(n_genes,1);  % An array for each gene.
% (i,j) is the combination which individual i has in clique j (in haploid case..).

for i = 1:n_genes
    if n_cliques(i)<1
        cliq_data{i} = zeros(ninds, 1);
        if gene_lengths(i)==2
            % linkage_model must be 'codon' to end up here..
            positions = component_mat(i, [1 2]);
            rows = data(:,positions);
            observations = (rows(:,1)-1) * noalle(positions(2)) + rows(:,2);
        else
            positions = component_mat(i,1);
            rows = data(:,positions);
            observations = rows;
        end
        cliq_data{i}(:,1) = observations;
    else
        cliq_data{i} = zeros(ninds, n_cliques(i));
        for j = 1:n_cliques(i)
            if strcmp(linkage_model,'codon')
                positions = component_mat(i,j:j+2);
                rows = data(:, positions);
                observations = (rows(:,1)-1) * prod(noalle(positions(2:3))) + ...
                    (rows(:,2)-1) * noalle(positions(3)) + rows(:,3);
            else
                positions = component_mat(i,j:j+1);
                rows = data(:,positions);
                observations = (rows(:,1)-1) * noalle(positions(2)) + rows(:,2);
            end
            cliq_data{i}(:,j) = observations;
        end
    end
end

cliq_counts = cell(n_genes,1);
% (i,j,k) is the count of combination i, in clique j, in population k.

npops = length(unique(PARTITION));
for i = 1:n_genes
    cliq_counts{i} = zeros(max_noalle(i), max(1,n_cliques(i)), npops);
    for j = 1:npops
        partition = repmat(PARTITION', [rowsFromInd 1]);
        partition = partition(:);  % Partition for rows in the data (instead of individuals).
        inds_now = find(partition==j);
        ninds_now = length(inds_now);
        data_now = cliq_data{i}(inds_now,:);

        for k = 1:max(n_cliques(i),1)
            apu = zeros(ninds_now, max_noalle(i));
            apu(sub2ind([ninds_now max_noalle(i)],...
                (1:ninds_now)', data_now(:,k)))=1;
            cliq_counts{i}(:, k, j) = (sum(apu,1))';
        end
    end
end

sep_data = cell(n_genes,1);
n_separators = n_cliques-1;
for i = 1:n_genes
    sep_data{i} = zeros(ninds, n_separators(i));
    for j = 1:n_separators(i)
        if strcmp(linkage_model, 'codon')
            positions = component_mat(i,j+1:j+2);
            rows = data(:, positions);
            observations = (rows(:,1)-1) * noalle(positions(2)) + rows(:,2);
        else
            positions = component_mat(i,j+1);
            rows = data(:,positions);
            observations = rows;
        end
        sep_data{i}(:,j) = observations;
    end
end


%------------------------------------------------------


function [cliq_freqs, sep_freqs] = simulateCliqFreqs(cliq_counts, noalle, component_mat, ...
    gene_lengths, linkage_model)

% cliq_freqs: cell-array. Element (i,j,k) in cell m is the frequence of
% combination i, in clique j of the m:th gene, in population k.

% sep_freqs: like cliq_freqs, but for the separators.

%------------------------------------------------------------------------

% cliq_freqs: cell-array, jossa on vastaavat dimensiot kuin cliq_counts:issa.
% solun m alkio (i,j,k) on geenin m, klikin j, kombinaation i, frekvenssi
% populaatiossa k.

% sep_freqs: cell-array, kuten cl_freqs, mutta separaattoreille.

%-------------------------------------------------------------------------

global PARTITION;

n_genes = length(cliq_counts);
if strcmp(linkage_model,'codon')
    n_cliques = gene_lengths-2;  % Number of cliques in each gene
else
    n_cliques = gene_lengths-1;  % Use linear model
end
npops = length(unique(PARTITION));

cliq_freqs = cell(n_genes,1);
sep_freqs = cell(n_genes,1);

for i=1:n_genes
    
    cliq_freqs{i} = zeros(size(cliq_counts{i}));
    
    positions = component_mat(i,1:gene_lengths(i));
    
    if n_cliques(i)<1
        % the gene is shorter than a normal clique.
        if gene_lengths(i)==2
            % linkage_model must be 'codon' to end up here..
            cliq_noalle = noalle(positions(1)) .* noalle(positions(2));
            sep_noalle = [];
        else
            % gene_lengths(i) == 1
            cliq_noalle = noalle(positions(1));
            sep_noalle = [];
        end
        sep_freqs{i} = [];
    else
        if strcmp(linkage_model, 'codon')
            cliq_noalle = noalle(positions(1:end-2)) .* noalle(positions(2:end-1)) .* ...
                noalle(positions(3:end));
            sep_noalle = noalle(positions(2:end-2)) .* noalle(positions(3:end-1));
        else
            cliq_noalle = noalle(positions(1:end-1)) .* noalle(positions(2:end));
            sep_noalle = noalle(positions(2:end-1));
        end
        sep_freqs{i} = zeros(max(sep_noalle), n_cliques(i)-1, npops);
    end
    
    % First clique:
    prior = (1 / cliq_noalle(1)) * ones(cliq_noalle(1),1);
    counts_now = repmat(prior, [1 1 npops]) + cliq_counts{i}(1:cliq_noalle(1),1,:);
    for k = 1:npops
        simul = randdir(counts_now(:,1,k), cliq_noalle(1));
        cliq_freqs{i}(1:cliq_noalle(1),1,k) = simul;
    end
    
    for j=2:n_cliques(i)
        % Obtain freqs for j-1:th separator by marginalization from j-1:th
        % clique, and draw values for the frequencies of the j:th clique:
        
        aux = cliq_freqs{i}(1:cliq_noalle(j-1), j-1, :); % Freqs of the previous clique
        aux = reshape(aux, [sep_noalle(j-1), noalle(positions(j-1)), npops]);
        
        % Freqs for separator by marginalization from the previous clique:
        sep_freqs{i}(1:sep_noalle(j-1),j-1,:) = sum(aux,2);
            
        prior = (1 / cliq_noalle(j)) * ones(cliq_noalle(j),1);
        counts_now = repmat(prior, [1 1 npops]) + cliq_counts{i}(1:cliq_noalle(j),j,:);
        for k = 1:npops
            % Simulate conditional frequencies:
            for m = 1:sep_noalle(j-1)
                if strcmp(linkage_model, 'codon')
                    values = (m-1)*noalle(positions(j+2))+1:m*noalle(positions(j+2));
                else
                    values = (m-1)*noalle(positions(j+1))+1:m*noalle(positions(j+1));
                end
                simul = randdir(counts_now(values,1,k), length(values));
                cliq_freqs{i}(values,j,k) = simul * sep_freqs{i}(m,j-1,k);  % MIETI TARKKAAN!
            end
        end
    end
end


%--------------------------------------------------------------------


function [ownCliqFreqs, ownSepFreqs] = computePersonalCliqueFreqs(...
    ind, cl_data, cl_freqs, sep_data, sep_freqs, rowsFromInd, ...
    gene_lengths, linkage_model)

% ownCliqFreqs is (npops * (n_cliques*rowsFromInd)) table, where each column
% contains the frequencies of the corresponding clique_combination, in 
% different populations.

% ownSepFreqs is (npops * (n_seps*rowsFromInd)) table, like ownCliqFreqs.

n_genes = length(gene_lengths);
if strcmp(linkage_model,'codon')
    n_cliques = gene_lengths-2;
    n_cliques = max([n_cliques ones(n_genes,1)], [], 2);  % for genes shorter than clique
else
    n_cliques = gene_lengths-1;  % Use linear model.
    n_cliques = max([n_cliques ones(n_genes,1)], [], 2);
end

total_n_cliques = sum(n_cliques);
npops = size(cl_freqs{1},3);

ownCliqFreqsXX = zeros(1, total_n_cliques*rowsFromInd, npops);

pointer = 1;
for i = 1:n_genes
    ind_data = cl_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd , :);
    for j = 1:n_cliques(i)   % MUUTA!
        for k = 1:rowsFromInd
            code = ind_data(k,j);
            ownCliqFreqsXX(1,pointer,:) = cl_freqs{i}(code,j,:);
            pointer = pointer+1;
        end
    end
end

ownCliqFreqs = (squeeze(ownCliqFreqsXX))';

n_separators = n_cliques-1;
total_n_separators = sum(n_separators);

ownSepFreqsXX = zeros(1, total_n_separators*rowsFromInd, npops);

pointer = 1;
for i = 1:n_genes
    ind_data = sep_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd , :);
    for j = 1:n_separators(i)
        for k = 1:rowsFromInd
            code = ind_data(k,j);
            ownSepFreqsXX(1,pointer,:) = sep_freqs{i}(code,j,:);
            pointer = pointer+1;
        end
    end
end

if (total_n_separators*rowsFromInd)==1
    ownSepFreqs = (squeeze(ownSepFreqsXX));
else
    ownSepFreqs = (squeeze(ownSepFreqsXX))';
end


%-------------------------------------------------------------------------


function exp_cliq_freqs = computeExpectedFreqs(cliq_counts, ...
    noalle, component_mat, gene_lengths, linkage_model)

% Returns the expected values for the clique and separator frequencies in
% different populations.

% exp_cliq_freqs: cell-array. Element (i,j,k) in cell m is the expected
% frequence of combination i, in clique j of the m:th gene, in 
% population k.

n_genes = length(gene_lengths);

if strcmp(linkage_model, 'codon')
    n_cliques = gene_lengths-2;
else
    n_cliques = gene_lengths-1;  % Linear model
end

npops = size(cliq_counts{1},3);
exp_cliq_freqs = cell(n_genes,1);

for i = 1:n_genes
    
    exp_cliq_freqs{i} = zeros(size(cliq_counts{i}));
    positions = component_mat(i,1:gene_lengths(i));
    
    if n_cliques(i)<1
        % the gene is shorter than a normal clique.
        if gene_lengths(i)==2
            % linkage_model must be 'codon' to end up here..
            cliq_noalle = noalle(positions(1)) .* noalle(positions(2));
        else
            % gene_lengths(i) == 1
            cliq_noalle = noalle(positions(1));
        end
    else
        if strcmp(linkage_model, 'codon')
            cliq_noalle = noalle(positions(1:end-2)) .* noalle(positions(2:end-1)) .* ...
                noalle(positions(3:end));
        else
            cliq_noalle = noalle(positions(1:end-1)) .* noalle(positions(2:end));
        end
    end
    
    for j = 1:max(1, n_cliques(i))
        prior = (1 / cliq_noalle(j)) * ones(cliq_noalle(j),1);
        counts_now = repmat(prior, [1 1 npops]) + cliq_counts{i}(1:cliq_noalle(j),j,:);
        exp_cliq_freqs{i}(1:cliq_noalle(j),j,:) = ...
            counts_now ./ repmat(sum(counts_now,1), [cliq_noalle(j) 1 1]);
    end
end


%----------------------------------------------------------


function [ref_cliq_data, ref_sep_data] = ...
    simulateLinkageIndividuals(n, rowsFromInd, exp_cliq_freqs, gene_lengths, ...
    noalle, component_mat, linkage_model)

% Simulates n individuals from each population using expected frequencies 
% for cliques and separators.

% ref_cliq_data: cell array, each cell corresponds to one gene. Elements
% ((i-1)*rowsFromInd+1:i*rowsFromInd, j) in cell k are the codes of the allele 
% combinations in the j:th clique in gene k, for individual i.

% ref_sep_data: like ref_cliq_data. i:th separator separates cliques i and i+1.

n_genes = length(gene_lengths);
if strcmp(linkage_model,'codon')
    n_cliques = gene_lengths-2;
else
    n_cliques = gene_lengths-1;  % Linear model
end
npops = size(exp_cliq_freqs{1},3);
ninds = n*npops;

ref_cliq_data = cell(n_genes,1);
ref_sep_data = cell(n_genes,1);

for i = 1:n_genes
    ref_cliq_data{i} = zeros(ninds*rowsFromInd, max(n_cliques(i),1));   % Added: rowsFromInd
    
    positions = component_mat(i,1:gene_lengths(i));
    
    if strcmp(linkage_model,'codon')
        sep_noalle = noalle(positions(2:end-2)) .* noalle(positions(3:end-1));
    else
        sep_noalle = noalle(positions(2:end-1));
    end
    ref_sep_data{i} = zeros(ninds*rowsFromInd, n_cliques(i)-1);  % Added: rowsFromInd
    
    for ind = 1:ninds
        pop = ceil(ind/n);
        
        % First clique:
        freqs = exp_cliq_freqs{i}(:,1,pop);
        freqs = repmat(freqs, [1 rowsFromInd]);
        codes = simulateCodes(freqs);
        ref_cliq_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd, 1) = codes;
        
        for j = 2:n_cliques(i)
            previous_cliq = ref_cliq_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd, j-1);
            
            % Value for j-1:th separator:
            sep_codes = rem(previous_cliq, sep_noalle(j-1));
            sep_codes(find(sep_codes==0)) = sep_noalle(j-1);
            ref_sep_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd, j-1) = sep_codes;
            
            % Value for j:th clique:
            if strcmp(linkage_model,'codon')
                freqs = zeros(noalle(positions(j+2)),rowsFromInd);
                values = zeros(noalle(positions(j+2)),rowsFromInd);
                for k = 1:rowsFromInd
                    values(:,k) = ((sep_codes(k)-1)*noalle(positions(j+2))+1 : sep_codes(k)*noalle(positions(j+2)))';
                    freqs(:,k) = exp_cliq_freqs{i}(values(:,k), j, pop);
                end
                freqs = freqs ./ repmat(sum(freqs,1), [noalle(positions(j+2)) 1]);
            else
                freqs = zeros(noalle(positions(j+1)),rowsFromInd);
                values = zeros(noalle(positions(j+1)),rowsFromInd);
                for k = 1:rowsFromInd
                    values(:,k) = ((sep_codes(k)-1)*noalle(positions(j+1))+1 : sep_codes(k)*noalle(positions(j+1)))';
                    freqs(:,k) = exp_cliq_freqs{i}(values(:,k), j, pop);
                end
                freqs = freqs ./ repmat(sum(freqs,1), [noalle(positions(j+1)) 1]);
            end
            codes = simulateCodes(freqs);
            codes = values(sub2ind(size(values),codes,(1:rowsFromInd)'));
            ref_cliq_data{i}((ind-1)*rowsFromInd+1:ind*rowsFromInd, j) = codes;
        end
    end
end

function codes = simulateCodes(freqs)
% Freqs is a table where each column is a distribution. The
% number of columns in freqs must be equal to rowsFromInd.
% A value is drawn from each distribution in different columns. The values
% are saved in codes, which is (rowsFromInd*1) table.

[nrows, rowsFromInd] = size(freqs);
codes = nrows+1-sum(cumsum(freqs)>repmat(rand(1,rowsFromInd),[nrows 1]),1);
codes = codes';