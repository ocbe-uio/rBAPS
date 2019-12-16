function semi_linkageMixture_speed(c_train, c_test)
% This function process adjusts the priors of the training data accoring to
% the test data. Based on the adjusted priors, the test data is clustered. 

% modified from linkageMixture_speed.m by Lu Cheng, 16.02.2010

% Update by Lu Cheng, 07.03.2011
% case of only 1 sample in the test data has been handled

% added by Lu Cheng, 11.03.2010
global SCRIPT_MODE;
global PARAMETERS;
if isempty(SCRIPT_MODE)
    SCRIPT_MODE = false;
end
% -----------------


%% compare the training data and test data, adjust priors

%1% Compare the training data and test data to adjust the prior
if ~all(all(c_train.component_mat == c_test.component_mat))
    disp('The gene lengths are different between the training data and the test data!');
    return;
end

flag = false; % whether the trained priors should be adjusted
n_loci = size(c_train.alleleCodes,2);

if c_train.rowsFromInd ~= c_test.rowsFromInd
    error('Inconsistant rows from each individual. Train: %d Test: %d. Quit! \n', ...
        c_train.rowsFromInd, c_test.rowsFromInd);
    return;
elseif c_train.rowsFromInd > 1
    error('Data must be haploid. Quit! rowsFromInd: %d.\n', c_train.rowsFromInd);
    return;
end

for i=1:n_loci
    if flag;   break;    end
    a = setdiff(c_test.alleleCodes(:,i),c_train.alleleCodes(:,i));
    a = a(a~=0);
    if ~isempty(a)
        flag = true;
        fprintf('New alleles are detected in the test data at loci: %d\n',i);
        fprintf('The processing time will be much longer than usual. \n');
    end
end

%2% reprocess the data for clustering
if flag
    
    % combing the training data and test data, adjust the priors
    combine_data = [c_train.data; c_test.data];
    n_train = size(c_train.data,1);
    n_samples = size(combine_data,1);
    combine_data(:,end) = (1:n_samples)';
    [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(combine_data);
    
    index = data(:,end);
    [data_clique, data_separator, noalle_clique, noalle_separator] = ...
        transform4(data, c_train.component_mat,'codon');
    data_clique = [data_clique index];
    data_separator = [data_separator index];
    
    % Count the data
    [counts_cq, nalleles_cq, prior_cq, adjprior_cq, genotypes_cq]...
        = allfreqsnew2(data_clique, double(noalle_clique));
    [counts_sp, nalleles_sp, prior_sp, adjprior_sp, genotypes_sp]...
        = allfreqsnew2(data_separator, double(noalle_separator));

    counts_cq = uint16(counts_cq);
    counts_sp = uint16(counts_sp);
        
    c_train.adjprior = adjprior;
    
    c_train.counts_cq = counts_cq(:,:,1:n_train);
    c_train.counts_sp = counts_sp(:,:,1:n_train);
    
    c_test.counts_cq = counts_cq(:,:,n_train+1:end);
    c_test.counts_sp = counts_sp(:,:,n_train+1:end);
    
    c_train.adjprior_cq = adjprior_cq;
    c_train.adjprior_sp = adjprior_sp;
    
    c_train.alleleCodes = alleleCodes;
    c_train.noalle = noalle;
        
    clear data rowsFromInd alleleCodes noalle adjprior priorTerm index 
    clear data_clique data_separator noalle_clique noalle_separator
    clear counts_cq counts_sp nalleles_cq nalleles_sp prior_cq prior_sp adjprior_cq adjprior_sp
    
else
    
    % KEY: adjust the test data to fit the configuration of the training data
    %      the 'codes_cq' and 'codes_sp' are directly translated from DNA sequence
    %      SEE 'i_encode_n.m' under the linkage folder
    
    num_cq  = size(c_train.counts_cq, 1);   num_sp  = size(c_train.counts_sp, 1); 
    n_loci_cq = size(c_train.counts_cq, 2); n_loci_sp = size(c_train.counts_sp, 2);
    n_inds = size(c_test.counts_cq, 3);
    
    counts_cq = zeros(num_cq, n_loci_cq, n_inds);
    counts_sp = zeros(num_sp, n_loci_sp, n_inds);
    
    % mapping the indexes of cliques and separators of the test data to the
    % indexes of the training data
    for k = 1:n_inds
        for j = 1:n_loci_cq
            [c, ia, ib] = intersect(c_test.codes_cq{j}, c_train.codes_cq{j},'rows');
            counts_cq(ib,j,k) = c_test.counts_cq(ia,j,k);
        end
        
        for j = 1:n_loci_sp
            [c, ia, ib] = intersect(c_test.codes_sp{j}, c_train.codes_sp{j},'rows');
            counts_sp(ib,j,k) = c_test.counts_sp(ia,j,k);
        end
    end

    c_test.counts_cq = counts_cq;
    c_test.counts_sp = counts_sp;

    clear c ia ib k j i;

end

%% cluster the test data

% case of only 1 sample in the test data, added by Lu Cheng, 07.03.2011
if size(c_test.data,1)~=1
    [Z,dist] = newGetDistances(c_test.data, c_test.rowsFromInd);
    c_test.Z = Z;
    c_test.dist = dist;
end

clear Z dist;

message = cat(2,'There are currently ',num2str(length(unique(c_train.cluster_labels))),' clusters in the training data, please input upper bounds of cluster numbers in the test data.');

if SCRIPT_MODE
    cluster_nums = str2num(PARAMETERS.cluster_num_upperbounds);
else
    cluster_nums = inputdlg(message);
    if isempty(cluster_nums) == 1
        return;
    else
        cluster_nums = str2num(cluster_nums{:});
    end    
end

% % Test purpose, Check the input data, there should be 1 allele in each loci (column)
% % Lu Cheng, 25.02.2010
% if ~all(all(squeeze(sum(c_train.counts_cq,1))))
%     disp('Missing cq value of some sample in counts_cq of the training data');
%     return;
% elseif ~all(all(squeeze(sum(c_train.counts_sp,1))))
%     disp('Missing sp value of some sample in counts_sp of the training data');
%     return;
% elseif ~all(all(squeeze(sum(c_test.counts_cq,1))))
%     disp('Missing cq value of some sample in counts_cq of the test data');
%     return;
% elseif ~all(all(squeeze(sum(c_test.counts_sp,1))))
%     disp('Missing sp value of some sample in counts_sp of the test data');
%     return;
% end

tic
semi_res = semi_linkageMix(c_train, c_test, cluster_nums);
toc

semi_res.popnames = c_test.popnames;

writeMixtureInfo(semi_res);

% save the results 
if SCRIPT_MODE
    save_results = PARAMETERS.save_results;
else
    save_results = questdlg('Do you wish to save the results?',...
           'Save Results','Yes','No','Yes');
end

if isequal(save_results,'Yes')
    if SCRIPT_MODE
        save(PARAMETERS.result_file,'semi_res','-v7.3');
    else
        [filename, pathname] = uiputfile('*.mat','Save the results as');
        if (sum(filename)==0) || (sum(pathname)==0)
            % do nothing
        else
            save(strcat(pathname,filename),'semi_res','-v7.3');
        end
    end   
end;
        
% -----------------------------------------------------------------------



%--------------------------------------------------------------------------
%% The next three functions are for computing the initial partition
% according to the distance between the individuals

function initial_partition=admixture_initialization(nclusters,Z)
T=cluster_own(Z,nclusters);
initial_partition=T;

function T = cluster_own(Z,nclust)
% true=logical(1);
% false=logical(0);

maxclust = nclust;
% Start of algorithm
m = size(Z,1)+1;
T = zeros(m,1);
% maximum number of clusters based on inconsistency
if m <= maxclust
    T = (1:m)';
elseif maxclust==1
    T = ones(m,1);
    
else
    clsnum = 1;
    for k = (m-maxclust+1):(m-1)
        i = Z(k,1); % left tree
        if i <= m % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
        i = Z(k,2); % right tree
        if i <= m  % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
    end
end

function T = clusternum(X, T, k, c)
m = size(X,1)+1;
while(~isempty(k))
    % Get the children of nodes at this level
    children = X(k,1:2);
    children = children(:);

    % Assign this node number to leaf children
    t = (children<=m);
    T(children(t)) = c;

    % Move to next level
    k = children(~t) - m;
end

%--------------------------------------------------------------------------

function Z = computeLinkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 || m ~= fix(m)
    error('The first input has to match the output of the PDIST function in size.');
end
if nargin == 1 % set default switch to be 'co'
    method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
% monotonic = 1;
Z = zeros(m-1,3); % allocate the output matrix.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;
for s = 1:(n-1)
    X = Y;
    [v, k] = min(X);
    i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
    j = k - (i-1)*(m-i/2)+i;
    Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
    I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
    U = [I1 I2 I3];
    I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
    J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

    switch method
        case 'si' %single linkage
            Y(I) = min(Y(I),Y(J));
        case 'av' % average linkage
            Y(I) = Y(I) + Y(J);
        case 'co' %complete linkage
            Y(I) = max(Y(I),Y(J));
        case 'ce' % centroid linkage
            K = N(R(i))+N(R(j));
            Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
        case 'wa'
            Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
                N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
    end
    J = [J i*(m-(i+1)/2)-m+j];
    Y(J) = []; % no need for the cluster information about j.

    % update m, N, R
    m = m-1;
    N(n+s) = N(R(i)) + N(R(j));
    R(i) = n+s;
    R(j:(n-1))=R((j+1):n);
end

%--------------------------------------------------------------------------

%-------------------------------------------------------------------------

function [Z, dist] = newGetDistances(data, rowsFromInd)

ninds = max(data(:,end));
nloci = size(data,2)-1;
riviLkm = nchoosek(double(ninds),2);

% empties = find(data<0);
% data(empties)=0;
data(logical(data<0)) = 0;
data = uint16(data);

pariTaulu = zeros(riviLkm,2);
aPointer=1;

for a=1:ninds-1
    pariTaulu(aPointer:aPointer+double(ninds-1-a),1) = ones(ninds-a,1,'uint16')*a;
    pariTaulu(aPointer:aPointer+double(ninds-1-a),2) = uint16((a+1:ninds)');
    aPointer = aPointer+double(ninds-a);
end

eka = pariTaulu(:,ones(1,rowsFromInd));
eka = eka * rowsFromInd;
miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
eka = eka - miinus;

toka = pariTaulu(:,ones(1,rowsFromInd)*2);
toka = toka * rowsFromInd;
toka = toka - miinus;

eka = uint16(eka);
toka = uint16(toka);

clear pariTaulu; clear miinus;

summa = uint16(zeros(riviLkm,1));
vertailuja = uint16(zeros(riviLkm,1));

x = zeros(size(eka));    x = uint16(x);
y = zeros(size(toka));   y = uint16(y);
% fprintf(1,'%%10');
for j=1:nloci;
    
    for k=1:rowsFromInd
        x(:,k) = data(eka(:,k),j);
        y(:,k) = data(toka(:,k),j);
    end

    for a=1:rowsFromInd
        for b=1:rowsFromInd
            vertailutNyt = uint16(x(:,a)>0 & y(:,b)>0);
            vertailuja = vertailuja + vertailutNyt;
            lisays = (x(:,a)~=y(:,b) & vertailutNyt);
            summa = summa + uint16(lisays);
        end
    end
    % fprintf(1,'\b\b');
    % fprintf(1,'%d',floor(10+80*j/nloci));
end

clear x;    clear y;   clear vertailutNyt;
clear eka; clear toka; clear data; clear lisays;
dist = zeros(length(vertailuja),1);
% nollat = find(vertailuja==0);
% dist(nollat) = 1;
dist(logical(vertailuja==0)) = 1;
muut = find(vertailuja>0);
dist(muut) = double(summa(muut))./double(vertailuja(muut));
clear summa; clear vertailuja; clear muut;

Z = computeLinkage(dist');
% fprintf(1,'\b\b');
% fprintf(1,'%d\n',100);
%--------------------------------------------------------------------------

function writeMixtureInfo(c)

outputFile = 'baps5_semi_output.txt';

% output the semi-supervised clustering results to the outputFile
% modified by Lu Cheng, 28.03.2010

ninds = length(c.PARTITION);
npops =  c.npops;
popnames = c.popnames;
logml = c.logml;
partition = c.PARTITION;
partitionSummary = c.partitionSummary;

if ~isempty(outputFile)
    fid = fopen(outputFile,'w');
else
    fid = -1;
    %diary('baps5_semi_output.baps'); % save in text anyway.
end

dispLine;
disp('RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:');
disp(['Number of clustered individuals: ' ownNum2Str(ninds)]);
disp(['Number of groups in optimal partition: ' ownNum2Str(npops)]);
disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
    fprintf(fid,'%10s\n', ['RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:']);
    fprintf(fid,'%20s\n', ['Number of clustered individuals: ' ownNum2Str(ninds)]);
    fprintf(fid,'%20s\n', ['Number of groups in optimal partition: ' ownNum2Str(npops)]);
    fprintf(fid,'%20s\n\n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
end

disp('Best Partition: ');
if (fid ~= -1)
    fprintf(fid,'%s \n','Best Partition: ');
end
for m=1:npops
    indsInM = find(partition==m);
    
    if isempty(indsInM)
        continue;
    end
    
    length_of_beginning = 11 + floor(log10(m));
    cluster_size = length(indsInM);
    
    text = ['Cluster ' num2str(m) ': {' char(popnames{indsInM(1)})];
    for k = 2:cluster_size
        text = [text ', ' char(popnames{indsInM(k)})];
    end;
    text = [text '}'];
    
    while length(text)>58
        %Take one line and display it.
        new_line = takeLine(text,58);
        text = text(length(new_line)+1:end);
        disp(new_line);
        if (fid ~= -1)
            fprintf(fid,'%s \n',new_line);
        end
        if length(text)>0
            text = [blanks(length_of_beginning) text];
        else
            text = [];
        end;
    end;
    
    if ~isempty(text)
        disp(text);
        if (fid ~= -1)
            fprintf(fid,'%s \n',text);
        end
    end;
end

names = true;

clusterProbTable = c.clusterProbTable;
if npops == 1
    clusterProbTable = [];
else
    disp('');
    disp('Posterior probability of assignment into clusters:');
    
    if (fid ~= -1)
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', 'Posterior probability of assignment into clusters: '); fprintf(fid, '\n');
    end

    text = sprintf('%10s','ind');
    for ii = 1:npops
        tmpstr = sprintf('%10s',num2str(ii));
        text = [text tmpstr];
    end
    
    disp(text);
    if (fid ~= -1)
        fprintf(fid, '%s \n', text);
    end
        
    for ii = 1:ninds
        text = sprintf('%10s',popnames{ii}{:});
        for jj = 1:npops
            tmpstr = sprintf('%10s',num2str(clusterProbTable(ii,jj),'%10.6f'));
            text = [text tmpstr];
        end
        
        if ii<100
            disp(text);
        elseif ii==101
            disp('.......................................');
            disp('..........see output file..............');
        end
        if (fid ~= -1)
            fprintf(fid, '%s \n', text);
        end   
        text = [];
    end
end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
    fprintf(fid, '%s \n\n', ' ');
    fprintf(fid, '%s \n', 'List of sizes of 10 best visited partitions and corresponding log(ml) values'); fprintf(fid, '\n');
end

partitionSummary = sortrows(partitionSummary,2);
partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
partitionSummary = partitionSummary(logical(partitionSummary(:,2)>-1e49),:);
if size(partitionSummary,1)>10
    vikaPartitio = 10;
else
    vikaPartitio = size(partitionSummary,1);
end
for part = 1:vikaPartitio
    line = [num2str(partitionSummary(part,1),'%20d') '    ' num2str(partitionSummary(part,2),'%20.6f')];
    disp(line);
    if (fid ~= -1)
        fprintf(fid, '%s \n', line);
    end
end

if (fid ~= -1)
    fclose(fid);
else
    diary off
end
    
%--------------------------------------------------------------
function newline = takeLine(description,width)
%Returns one line from the description: line ends to the first
%space after width:th mark.
% newLine = description(1:width);
n = width+1;
while ~isspace(description(n)) && n<length(description)
    n = n+1;
end;
newline = description(1:n);


function dispLine
disp('---------------------------------------------------');

function num2 = omaRound(num)
% Pyï¿½ristï¿½ï¿½ luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------

%-------------------------------------------------------------------------

function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
    handleData(raw_data)
% Alkuperï¿½isen datan viimeinen sarake kertoo, miltï¿?yksilï¿½ltï¿?
% kyseinen rivi on perï¿½isin. Funktio tutkii ensin, ettï¿?montako
% riviï¿?maksimissaan on perï¿½isin yhdeltï¿?yksilï¿½ltï¿? jolloin saadaan
% tietï¿½ï¿½ onko kyseessï¿?haploidi, diploidi jne... Tï¿½mï¿½n jï¿½lkeen funktio
% lisï¿½ï¿½ tyhjiï¿?rivejï¿?niille yksilï¿½ille, joilta on perï¿½isin vï¿½hemmï¿½n
% rivejï¿?kuin maksimimï¿½ï¿½rï¿?
%   Mikï¿½li jonkin alleelin koodi on =0, funktio muuttaa tï¿½mï¿½n alleelin
% koodi pienimmï¿½ksi koodiksi, joka isompi kuin mikï¿½ï¿½n kï¿½ytï¿½ssï¿?oleva koodi.
% Tï¿½mï¿½n jï¿½lkeen funktio muuttaa alleelikoodit siten, ettï¿?yhden lokuksen j
% koodit saavat arvoja vï¿½lillï¿?1,...,noalle(j).

data = raw_data;
nloci=size(raw_data,2)-1;

dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
    isoinAlleeli = max(max(dataApu));
    dataApu(nollat) = isoinAlleeli+1;
    data(:,1:nloci) = dataApu;
end
% dataApu = []; 
% nollat = []; 
% isoinAlleeli = [];

noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
   %alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(logical(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
end
alleleCodes = zeros(max(noalle),nloci);
for i=1:nloci
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
end

for loc = 1:nloci
    for all = 1:noalle(loc)
        % data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
        data(logical(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;

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
for ind=lessThanMax'    %Kï¿½y lï¿½pi ne yksilï¿½t, joilta puuttuu rivejï¿?
    miss = maxRowsFromInd-rowsFromInd(ind);  % Tï¿½ltï¿?yksilï¿½ltï¿?puuttuvien lkm.
    for j=1:miss
        rowToBeAdded = emptyRow;
        rowToBeAdded(end) = ind;
        data(nrows+pointer, :) = rowToBeAdded;
        pointer = pointer+1;
    end
end
data = sortrows(data, ncols);   % Sorttaa yksilï¿½iden mukaisesti
newData = data;
rowsFromInd = maxRowsFromInd;

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
end

%--------------------------------------------------------------------------
