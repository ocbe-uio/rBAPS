function c = preprocessXLS(xlsfile,varargin)
% This function preprocesses the input xlsfile
% File structure:  first line - title
%                  else - first column, name of individuals
%                       - second to end column, sequences of the given genes

% Lu Cheng, 16.02.2010

%% check file names

file_suff = '.xls'; % endings of the input file
if ~(exist(xlsfile,'file')==2)
    fprintf('Input file %s does not exists, quit!\n',xlsfile);
    return;
end

if ~strcmp(xlsfile(end-length(file_suff)+1:end),file_suff)
    fprintf('Input file %s does not end with %s, quit!\n',xlsfile,file_suff);
    return;
end

%% process the xls file

% Here we assume there is no missing values, if so, the missing values are
% indicated by 0
[data, component_mat, popnames] = processxls(xlsfile);

% missing data 0 is transformed to |alphabet|+1, in the [ACGT] case, '-' is 5
[data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);

c.data = data; c.rowsFromInd = rowsFromInd;
c.alleleCodes = alleleCodes; c.noalle=noalle;
c.adjprior = adjprior;
% c.priorTerm = c.priorTerm;

c.component_mat = component_mat;
c.popnames = popnames;

%% count the cliques and separators

index = data(:,end);

if isempty(varargin)
    [data_clique, data_separator, noalle_clique, noalle_separator, codes_cq, codes_sp, info_cq_loci, info_sp_loci] = ...
        transform5(data, component_mat);
else
    c_train = varargin{1};
    
    if ~all(all(c_train.component_mat == component_mat))
        disp('The gene lengths are different between the training data and the test data!');
        return;
    end
    
    [data_clique, data_separator, noalle_clique, noalle_separator, codes_cq, codes_sp, info_cq_loci, info_sp_loci] = ...
        transform5(data, component_mat, c_train.info_cq_loci,c_train.info_sp_loci);
end
data_clique = [data_clique index];
data_separator = [data_separator index];

% Count the data, note that the order of the alphabets keeps the same
[counts_cq, nalleles_cq, prior_cq, adjprior_cq, genotypes_cq] ...
    = allfreqsnew2(data_clique, double(noalle_clique));
[counts_sp, nalleles_sp, prior_sp, adjprior_sp, genotypes_sp] ...
    = allfreqsnew2(data_separator, double(noalle_separator));

clear prior_cq prior_sp nalleles_cq nalleles_sp genotypes_cq genotypes_sp;

counts_cq = uint16(counts_cq);
counts_sp = uint16(counts_sp);

c.counts_cq = counts_cq;
c.counts_sp = counts_sp;

c.adjprior_cq = adjprior_cq;
c.adjprior_sp = adjprior_sp;

c.codes_cq = codes_cq;
c.codes_sp = codes_sp;

c.info_cq_loci = info_cq_loci;
c.info_sp_loci = info_sp_loci;

%--------------------------------------------------------------------------
function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(raw_data)
% Alkuperäisen datan viimeinen sarake kertoo, miltä yksilöltä
% kyseinen rivi on peräisin. Funktio tutkii ensin, että montako
% riviä maksimissaan on peräisin yhdeltä yksilöltä, jolloin saadaan
% tietää onko kyseessä haploidi, diploidi jne... Tämän jälkeen funktio
% lisää tyhjiä rivejä niille yksilöille, joilta on peräisin vähemmän
% rivejä kuin maksimimäärä.
%   Mikäli jonkin alleelin koodi on =0, funktio muuttaa tämän alleelin
% koodi pienimmäksi koodiksi, joka isompi kuin mikään käytössä oleva koodi.
% Tämän jälkeen funktio muuttaa alleelikoodit siten, että yhden lokuksen j
% koodit saavat arvoja välillä 1,...,noalle(j).

% English Comments added 
% Small modifications have been added
% Lu Cheng, 17.02.2010

% Last column are the indexes of the samples, the raw_data is supposed to
% be unit16 type, 0 indicates missing value
data = raw_data;
nloci=size(raw_data,2)-1;

% Replace missing value with the |alphabet|+1, thus 0 is replaced by 5 for
% DNA dataset
dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
    isoinAlleeli = max(max(dataApu));
    dataApu(nollat) = isoinAlleeli+1;
    data(:,1:nloci) = dataApu;
end

% stores all different alleles at each loci, construct the allle codes matrix
noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(logical(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
end
alleleCodes = zeros(max(noalle),nloci);
for i=1:nloci
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
end

%-----------------modified by Lu Cheng 17.02.2010--------------------------%
% NOTE: Here we do not want to change the alpahbets, thus the following
% lines are commented

% replace the index of an allele to replace the allele
% for loc = 1:nloci
%     for all = 1:noalle(loc)
%         data(logical(data(:,loc)==alleleCodes(all,loc)), loc)=all;
%     end;
% end;
%-----------------modified end.....----------------------------------------%

% handle diploid situation
nind = max(data(:,end));
nrows = size(data,1);
ncols = size(data,2);
rowsFromInd = zeros(nind,1);
for i=1:nind
    rowsFromInd(i) = length(find(data(:,end)==i));
end
maxRowsFromInd = max(rowsFromInd);
a = -999;
emptyRow = repmat(a, 1, ncols);
lessThanMax = find(rowsFromInd < maxRowsFromInd);
missingRows = maxRowsFromInd*nind - nrows;
data = [data; zeros(missingRows, ncols)];
pointer = 1;
for ind=lessThanMax'    %Käy läpi ne yksilöt, joilta puuttuu rivejä
    miss = maxRowsFromInd-rowsFromInd(ind);  % Tältä yksilöltä puuttuvien lkm.
    for j=1:miss
        rowToBeAdded = emptyRow;
        rowToBeAdded(end) = ind;
        data(nrows+pointer, :) = rowToBeAdded;
        pointer = pointer+1;
    end
end
data = sortrows(data, ncols);   % Sorttaa yksilöiden mukaisesti
newData = data;
rowsFromInd = maxRowsFromInd;

% calculate the prior for each loci, priorTerm is a constant term in the
% formula, which is precalclulateed for speeding up the program
adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
end