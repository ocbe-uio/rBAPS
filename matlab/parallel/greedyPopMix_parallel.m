function greedyPopMix_parallel(options)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
clearGlobalVars;

% LASKENNAN ALKUARVOJEN MƒƒRITTƒMINEN
outp = [options.outputMat '.txt'];
inp = options.dataFile;

if isequal(options.dataType,'numeric')  %Raakadata
    data = load(options.dataFile);
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    if (ninds==0) 
        disp('*** ERROR: Incorrect Data-file.');
        return;    
    end
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    rowsFromInd = 0;  %Ei tiedet?
    
    if ~isempty(options.groupname)
        popnames = initPopNames(options.groupname);
        if (size(popnames,1)~=ninds)
            disp('*** ERROR: Incorrect name-file.');
            popnames = [];
        end
    else
        popnames = [];
    end
            
elseif isequal(options.dataType,'genepop')
    kunnossa = testaaGenePopData(options.dataFile);
    if kunnossa==0
        return
    end

    [data, popnames]=lueGenePopDataPop(options.dataType);
    [data, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(data);
    rowsFromInd = 2;  %Tiedet‰‰n GenePop:in tapauksessa.
       
end 

if ~isequal(options.dataType, 'matlab')
    a_data = data(:,1:end-1);

    npops = size(rows,1);
    PARTITION = 1:npops';  %Jokainen "yksil? eli populaatio on oma ryhm‰ns?
    [sumcounts, counts, logml] = ...
        initialPopCounts(a_data, npops, rows, noalle, adjprior);
    COUNTS = counts; SUMCOUNTS = sumcounts;
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);

    clear('counts', 'sumcounts','pathname','filename','vast2',...
        'vast3','vast4');
    [Z,dist] = getPopDistancesByKL(adjprior);  %Saadaan COUNTS:in avulla.
    
%     save_preproc = questdlg('Do you wish to save pre-processed data?',...
%         'Save pre-processed data?',...
%         'Yes','No','Yes');
%     if isequal(save_preproc,'Yes');
%         waitALittle;
%         [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
%         kokonimi = [pathname filename];
%         c.data = data; c.rows = rows; c.alleleCodes = alleleCodes;
%         c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
%         c.dist = dist; c.Z = Z; c.popnames = popnames; c.rowsFromInd = rowsFromInd;
%         c.npops = npops;  c.logml = logml;
%         save(kokonimi,'c');
%         clear c;
%     end;
end
        
if isequal(options.dataType, 'matlab')
    struct_array = load(options.dataFile);
    if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        if ~isfield(c,'rows')
            disp('Incorrect file format');
            return
        end
    elseif isfield(struct_array,'rows')  %Mideva versio
        c = struct_array;
    else
        disp('*** ERROR: Incorrect file format');
        return;
    end
    data = double(c.data); rows = c.rows; alleleCodes = c.alleleCodes;
    noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
    dist = c.dist; Z = c.Z; popnames = c.popnames; rowsFromInd = c.rowsFromInd;
    clear c;
end

if strcmp(options.fixedK, 'yes')
    fixedK = 1;
else
    fixedK = 0;
end

npopstext = [];
npopstextExtra = options.initialK;
if length(npopstextExtra)>=255
    npopstextExtra = npopstextExtra(1:255);
    npopstext = [npopstext ' ' npopstextExtra];
    teksti = 'The input field length limit (255 characters) was reached. Input more values: ';
else
    % -----------------------------------------------------
    % Set the limit of the input value.
    % Modified by Jing Tang, 30.12.2005
    if max(npopstextExtra) > size(data,1)
        error('Initial K larger than the sample size are not accepted. ');
    else
        npopstext = [npopstext ' ' num2str(npopstextExtra)];
    end
end

clear teksti;
if isempty(npopstext) || length(npopstext)==1
    return
else
    npopsTaulu = str2num(npopstext);
    ykkoset = find(npopsTaulu==1);
    npopsTaulu(ykkoset) = [];   % Mik‰li ykkˆsi?annettu yl‰rajaksi, ne poistetaan.
    if isempty(npopsTaulu)
        return
    end
    clear ykkoset;
end


c.data=data; c.rows = rows; c.alleleCodes = alleleCodes;
c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
c.dist=dist; c.Z=Z; c.rowsFromInd = rowsFromInd;

if fixedK
    % Only the first value of npopsTaulu is used
    npops = npopsTaulu(1);
    nruns = length(npopsTaulu);
    [logml, npops, partitionSummary]=indMix_fixK(c,npops,nruns,1);
else
    [logml, npops, partitionSummary]=indMix(c,npopsTaulu,1);   
end

if logml==1
    return
end

data = data(:,1:end-1);

viewPopMixPartition(PARTITION, rows, popnames);
%npops = poistaTyhjatPopulaatiot(npops);
%POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);

% h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
% h0 = findobj('Tag','filename2_text');
% outp = get(h0,'String');
changesInLogml = writeMixtureInfoPop(logml, rows, data, adjprior, priorTerm, ...
    outp,inp,partitionSummary, popnames, fixedK);
    
if exist('baps4_output.baps','file')
    copyfile('baps4_output.baps',outp)
    delete('baps4_output.baps')
end

if rowsFromInd==0
    %K‰ytettiin BAPS-formaattia, eik?rowsFromInd ole tunnettu.
    [popnames, rowsFromInd] = findOutRowsFromInd(popnames, rows);
end

groupPartition = PARTITION;

fiksaaPartitioYksiloTasolle(rows, rowsFromInd);
% The logml is saved for parallel computing
c.logml = logml;
c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
c.alleleCodes = alleleCodes; c.adjprior = adjprior;
c.rowsFromInd = rowsFromInd; c.popnames = popnames;
c.data = data; c.npops = npops; c.noalle = noalle;
c.mixtureType = 'popMix'; c.groupPartition = groupPartition;
c.rows = rows; c.changesInLogml = changesInLogml;
fprintf(1,'Saving the result...')
try
%     save(options.outputMat, 'c');
    save(options.outputMat, 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
    fprintf(1,'Finished.\n');
catch
    display('*** ERROR in saving the result.');
end


% -------------------------------------------------------------------------
% - Subfunctions
% -------------------------------------------------------------------------

function [newData, rows, alleleCodes, noalle, adjprior, priorTerm] = handlePopData(raw_data)
% Alkuper‰isen datan viimeinen sarake kertoo, milt?yksilˆlt?
% kyseinen rivi on per‰isin. Funktio muuttaa alleelikoodit 
% siten, ett?yhden lokuksen j koodit saavat arvoja 
% v‰lill?1,...,noalle(j). Ennen t‰t?muutosta alleeli, jonka
% koodi on nolla muutetaan.


data = raw_data;
nloci=size(raw_data,2)-1;

dataApu = data(:,1:nloci);
nollat = find(dataApu==0);
if ~isempty(nollat)
   isoinAlleeli = max(max(dataApu));
   dataApu(nollat) = isoinAlleeli+1;
   data(:,1:nloci) = dataApu;
end
dataApu = []; nollat = []; isoinAlleeli = [];

noalle=zeros(1,nloci);
alleelitLokuksessa = cell(nloci,1);
for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
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
        data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
end;

nind = max(data(:,end));
%rows = cell(nind,1);
rows = zeros(nind,2);
for i=1:nind
    rivit = find(data(:,end)==i)';
    rows(i,1) = min(rivit);
    rows(i,2) = max(rivit);
end
newData = data;

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j)); 
end


%----------------------------------------------------------------


function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];


%--------------------------------------------------------------------

function [Z,distances] = getPopDistancesByKL(adjprior)
% Laskee populaatioille et‰isyydet
% k‰ytt‰en KL-divergenssi?
global COUNTS;
maxnoalle = size(COUNTS,1);
nloci = size(COUNTS,2);
npops = size(COUNTS,3);
distances = zeros(nchoosek(npops,2),1);

d = zeros(maxnoalle, nloci, npops);
prior = adjprior;
prior(find(prior==1))=0;
nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
prior(1,nollia)=1;
for pop1 = 1:npops
    d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
    %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
end
pointer = 1;
for pop1 = 1:npops-1
    for pop2 = pop1+1:npops
        dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
        div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
        div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
        div = (div12+div21)/2;
        distances(pointer) = div;
        pointer = pointer+1;    
    end
end
Z=linkage(distances');

%--------------------------------------------------------------------------


function Z = linkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 | m ~= fix(m)
  error('The first input has to match the output of the PDIST function in size.');   
end
if nargin == 1 % set default switch to be 'co' 
   method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
monotonic = 1;
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


%-----------------------------------------------------------------------


function [sumcounts, counts, logml] = ...
    initialPopCounts(data, npops, rows, noalle, adjprior)

nloci=size(data,2);
counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);

for i=1:npops
    for j=1:nloci
        i_rivit = rows(i,1):rows(i,2);
        havainnotLokuksessa = find(data(i_rivit,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(i_rivit,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

logml = laskeLoggis(counts, sumcounts, adjprior);


%-----------------------------------------------------------------------


function loggis = laskeLoggis(counts, sumcounts, adjprior)
npops = size(counts,3);

logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior))) - ...
    sum(sum(gammaln(1+sumcounts)));
loggis = logml2;


%--------------------------------------------------------------------


function kunnossa = testaaGenePopData(tiedostonNimi)
% kunnossa == 0, jos data ei ole kelvollinen genePop data.
% Muussa tapauksessa kunnossa == 1.

kunnossa = 0;
fid = fopen(tiedostonNimi);
line1 = fgetl(fid);  %ensimm‰inen rivi
line2 = fgetl(fid);  %toinen rivi
line3 = fgetl(fid);  %kolmas

if (isequal(line1,-1) | isequal(line2,-1) | isequal(line3,-1))
    disp('Incorrect file format'); fclose(fid);
    return 
end
if (testaaPop(line1)==1 | testaaPop(line2)==1)
    disp('Incorrect file format'); fclose(fid);
    return
end
if testaaPop(line3)==1
    %2 rivi t‰llˆin lokusrivi
    nloci = rivinSisaltamienMjonojenLkm(line2);
    line4 = fgetl(fid);
    if isequal(line4,-1)
        disp('Incorrect file format'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin nelj?t‰ytyy sis‰lt‰‰ pilkku.
        disp('Incorrect file format'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedet‰‰n, ett?pys‰htyy
        pointer = pointer+1;
    end
    line4 = line4(pointer+1:end);  %pilkun j‰lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format'); fclose(fid);
        return
    end
else
    line = fgetl(fid);
    lineNumb = 4;
    while (testaaPop(line)~=1 & ~isequal(line,-1))
        line = fgetl(fid);
        lineNumb = lineNumb+1;    
    end
    if isequal(line,-1)
        disp('Incorrect file format'); fclose(fid);
        return
    end
    nloci = lineNumb-2;
    line4 = fgetl(fid);  %Eka rivi pop sanan j‰lkeen
    if isequal(line4,-1)
        disp('Incorrect file format'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin t‰ytyy sis‰lt‰‰ pilkku.
        disp('Incorrect file format'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedet‰‰n, ett?pys‰htyy.
        pointer = pointer+1;
    end
 
    line4 = line4(pointer+1:end);  %pilkun j‰lkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format'); fclose(fid);
        return
    end
end
kunnossa = 1;
fclose(fid);

%--------------------------------------------------------------------


function [data, popnames] = lueGenePopDataPop(tiedostonNimi)
% Data annetaan muodossa, jossa viimeinen sarake kertoo ryhm‰n.
% popnames on kuten ennenkin.

fid = fopen(tiedostonNimi);
line = fgetl(fid);  %ensimm‰inen rivi
line = fgetl(fid);  %toinen rivi
count = rivinSisaltamienMjonojenLkm(line);

line = fgetl(fid);
lokusRiveja = 1;
while (testaaPop(line)==0)
    lokusRiveja = lokusRiveja+1;
    line = fgetl(fid);
end

if lokusRiveja>1
    nloci = lokusRiveja;
else
    nloci = count;
end

popnames = cell(10,2);
data = zeros(100, nloci+1);
nimienLkm=0;
ninds=0;
poimiNimi=1;
digitFormat = -1;
while line ~= -1
    line = fgetl(fid);
    
    if poimiNimi==1  
        %Edellinen rivi oli 'pop'
        nimienLkm = nimienLkm+1;
        ninds = ninds+1;
        if nimienLkm>size(popnames,1);
            popnames = [popnames; cell(10,2)];
        end
        nimi = lueNimi(line);
        if digitFormat == -1
            digitFormat = selvitaDigitFormat(line);
            divider = 10^digitFormat;
        end
        popnames{nimienLkm, 1} = {nimi};   %N‰in se on greedyMix:iss‰kin?!?
        popnames{nimienLkm, 2} = ninds;
        poimiNimi=0;
        
        data = addAlleles(data, ninds, line, divider);
        
    elseif testaaPop(line)
        poimiNimi = 1;
        
    elseif line ~= -1
        ninds = ninds+1;
        data = addAlleles(data, ninds, line, divider);
    end
end

fclose(fid);
data = data(1:ninds*2,:);
popnames = popnames(1:nimienLkm,:);
npops = size(popnames,1);
ind = 1;
for pop = 1:npops
    if pop<npops
        while ind<popnames{pop+1,2}
            data([ind*2-1 ind*2],end) = pop;
            ind = ind+1;    
        end
    else
        while ind<=ninds
            data([ind*2-1 ind*2],end) = pop;
            ind = ind+1;    
        end
    end
end
        
%-------------------------------------------------------

function nimi = lueNimi(line)
%Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
n = 1;
merkki = line(n);
nimi = '';
while ~isequal(merkki,',')
    nimi = [nimi merkki];
    n = n+1;
    merkki = line(n);
end

%-------------------------------------------------------

function df = selvitaDigitFormat(line)
% line on ensimm‰inen pop-sanan j‰lkeinen rivi
% Genepop-formaatissa olevasta datasta. funktio selvitt‰‰
% rivin muodon perusteella, ovatko datan alleelit annettu
% 2 vai 3 numeron avulla.

n = 1;
merkki = line(n);
while ~isequal(merkki,',')
    n = n+1;
    merkki = line(n);
end

while ~any(merkki == '0123456789');
    n = n+1;
    merkki = line(n);
end
numeroja = 0;
while any(merkki == '0123456789');
    numeroja = numeroja+1;
    n = n+1;
    merkki = line(n);
end

df = numeroja/2;


%------------------------------------------------------


function count = rivinSisaltamienMjonojenLkm(line)
% Palauttaa line:n sis‰lt‰mien mjonojen lukum‰‰r‰n.
% Mjonojen v‰liss?t‰ytyy olla v‰lilyˆnti.
count = 0;
pit = length(line);
tila = 0;    %0, jos odotetaan v‰lilyˆntej? 1 jos odotetaan muita merkkej?
for i=1:pit
    merkki = line(i);
    if (isspace(merkki) & tila==0) 
        %Ei tehd?mit‰‰n.
    elseif (isspace(merkki) & tila==1)
        tila = 0;
    elseif (~isspace(merkki) & tila==0)
        tila = 1;
        count = count+1;
    elseif (~isspace(merkki) & tila==1)
        %Ei tehd?mit‰‰n
    end
end

%-------------------------------------------------------

function pal = testaaPop(rivi)
% pal=1, mik‰li rivi alkaa jollain seuraavista
% kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
% tapauksissa pal=0.

if length(rivi)<3 
    pal = 0;
    return
end
if (all(rivi(1:3)=='Pop') | ...
    all(rivi(1:3)=='pop') | ...
    all(rivi(1:3)=='POP'))
    pal = 1;    
    return
else 
    pal = 0;
    return
end


%--------------------------------------------------------


function data = addAlleles(data, ind, line, divider)
% Lisaa BAPS-formaatissa olevaan datataulukkoon
% yksilˆ‰ ind vastaavat rivit. Yksilˆn alleelit
% luetaan genepop-formaatissa olevasta rivist?
% line. Jos data on 3 digit formaatissa on divider=1000.
% Jos data on 2 digit formaatissa on divider=100.

nloci = size(data,2)-1;
if size(data,1) < 2*ind
    data = [data; zeros(100,nloci+1)];
end

k=1;
merkki=line(k);
while ~isequal(merkki,',')
   k=k+1;
   merkki=line(k);
end
line = line(k+1:end);
clear k; clear merkki;

alleeliTaulu = sscanf(line,'%d');

if length(alleeliTaulu)~=nloci
    disp('Incorrect data format.');
end

for j=1:nloci
    ekaAlleeli = floor(alleeliTaulu(j)/divider);
    if ekaAlleeli==0 ekaAlleeli=-999; end;
    tokaAlleeli = rem(alleeliTaulu(j),divider);
    if tokaAlleeli==0 tokaAlleeli=-999; end
    
    data(2*ind-1,j) = ekaAlleeli;
    data(2*ind,j) = tokaAlleeli;
end

data(2*ind-1,end) = ind;
data(2*ind,end) = ind;

%------------------------------------------------------------------------------------


function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on m‰‰ritelty pops-muuttujalla.

global COUNTS;
global SUMCOUNTS;
x = size(COUNTS,1);
y = size(COUNTS,2);
z = length(pops);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;


%--------------------------------------------------------------------------


function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, globalRows, data, adjprior, priorTerm)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li yksil?ind siirret‰‰n koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lis‰tt‰v?
% COUNTS:in siivuun i2, mik‰li muutos toteutetaan.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
muutokset = zeros(npops,1);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);

rows = globalRows(ind,1):globalRows(ind,2);
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%----------------------------------------------------------------------


function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukum‰‰r‰t (vastaavasti kuin COUNTS:issa), jotka ovat data:n 
% riveill?rows. rows pit‰‰ olla vaakavektori.

diffInCounts = zeros(max_noalle, nloci);
for i=rows
    row = data(i,:);
    notEmpty = find(row>=0);
    
    if length(notEmpty)>0
        diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
            diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
    end
end    

%------------------------------------------------------------------------


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, diffInCounts, ...
    adjprior, priorTerm)
% Suorittaa globaalien muuttujien muutokset, kun yksil?ind
% on siirret‰‰n koriin i2.

global PARTITION; 
global COUNTS; 
global SUMCOUNTS;
global POP_LOGML;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%--------------------------------------------------------------------------
%--

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2( ...
    i1, globalRows, data, adjprior, priorTerm);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li korin i1 kaikki yksilˆt siirret‰‰n
% koriin i. 

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
muutokset = zeros(npops,1);

i1_logml = POP_LOGML(i1);

inds = find(PARTITION==i1);
ninds = length(inds);

if ninds==0
    diffInCounts = zeros(size(COUNTS,1), size(COUNTS,2));
    return;
end

rows = [];
for i = 1:ninds
    ind = inds(i);
    lisa = globalRows(ind,1):globalRows(ind,2);
    rows = [rows; lisa'];
    %rows = [rows; globalRows{ind}'];
end

diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%---------------------------------------------------------------------------------


function updateGlobalVariables2( ...
    i1, i2, diffInCounts, adjprior, priorTerm);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilˆt siirret‰‰n koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;

inds = find(PARTITION==i1);
PARTITION(inds) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2, adjprior, priorTerm);


%--------------------------------------------------------------------------
%----

function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
    data, adjprior, priorTerm, i1)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mik?olisi muutos logml:ss? jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirret‰‰n koriin j.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
npops2 = length(unique(T2));
muutokset = zeros(npops2, npops);

i1_logml = POP_LOGML(i1);
for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
        rows = [];
        for i = 1:ninds
            ind = inds(i);
            lisa = globalRows(ind,1):globalRows(ind,2);
            rows = [rows; lisa'];
            %rows = [rows; globalRows{ind}'];
        end
        diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
        diffInSumCounts = sum(diffInCounts);

        COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
        new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
        COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';
    
        COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
        new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
        COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

        muutokset(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end    
end

%------------------------------------------------------------------------------------

function muutokset = laskeMuutokset5(inds, globalRows, data, adjprior, ...
    priorTerm, i1, i2)

% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mik‰li yksil?i vaihtaisi koria i1:n ja i2:n v‰lill?
    
global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;

ninds = length(inds);
muutokset = zeros(ninds,1);

i1_logml = POP_LOGML(i1);
i2_logml = POP_LOGML(i2);

for i = 1:ninds
    ind = inds(i);
    if PARTITION(ind)==i1
        pop1 = i1;  %mist?
        pop2 = i2;  %mihin
    else
        pop1 = i2;
        pop2 = i1;
    end
    rows = globalRows(ind,1):globalRows(ind,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)-diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)-diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)+diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)+diffInSumCounts;
    
    new_logmls = computePopulationLogml([i1 i2], adjprior, priorTerm);
    muutokset(i) = sum(new_logmls);
    
    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)+diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)+diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)-diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)-diffInSumCounts; 
end

muutokset = muutokset - i1_logml - i2_logml;

%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, ...
    adjprior, priorTerm, i2);
% Suorittaa globaalien muuttujien p‰ivitykset, kun yksilˆt 'muuttuvat'
% siirret‰‰n koriin i2. Ennen siirtoa yksilˆiden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


%----------------------------------------------------------------------------


function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sis‰lt‰‰ yksilˆiden inds2
% v‰liset et‰isyydet. ninds=kaikkien yksilˆiden lukum‰‰r?

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
rivi = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(rivi, 1) = inds2(i);
        apu(rivi, 2) = inds2(j);
        rivi = rivi+1;
    end
end
apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = dist(apu);


%-----------------------------------------------------------------------------------


function npops = poistaTyhjatPopulaatiot(npops)
% Poistaa tyhjentyneet populaatiot COUNTS:ista ja 
% SUMCOUNTS:ista. P‰ivitt‰‰ npops:in ja PARTITION:in.

global COUNTS;
global SUMCOUNTS;
global PARTITION;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);

for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
end
npops = length(notEmpty);


%-----------------------------------------------------------------------------------


function popnames = initPopNames(nameFile)

fid = fopen(nameFile);
if fid == -1
    %File didn't exist
    msgbox('Loading of the population names was unsuccessful', ...
        'Error', 'error');
    return;
end;
line = fgetl(fid);
counter = 1;
while (line ~= -1) & ~isempty(line)
    names{counter} = line;
    line = fgetl(fid);
    counter = counter + 1;
end;
fclose(fid);

popnames = cell(length(names), 2);
for i = 1:length(names)
    popnames{i,1} = names(i);
    popnames{i,2} = 0;
end


%-------------------------------------------------------------------------


function [popnames2, rowsFromInd] = findOutRowsFromInd(popnames, rows)

ploidisuus = questdlg('Specify the type of individuals in the data: ',...
    'Individual type?', 'Haploid', 'Diploid', 'Tetraploid', ...
    'Diploid');

switch ploidisuus    
case 'Haploid'
    rowsFromInd = 1;
case 'Diploid'
    rowsFromInd = 2;
case 'Tetraploid'
    rowsFromInd = 4;
end

if ~isempty(popnames)
    for i = 1:size(rows,1)
        popnames2{i,1} = popnames{i,1};
        rivi = rows(i,1):rows(i,2);
        popnames2{i,2} = (rivi(rowsFromInd))/rowsFromInd;    
    end
else 
    popnames2 = [];
end

%------------------------------------------------------------------

function fiksaaPartitioYksiloTasolle(rows, rowsFromInd)

global PARTITION;
totalRows = 0;
for ind = 1:size(rows,1)
    totalRows = totalRows + (rows(ind,2)-rows(ind,1)+1);
end
partitio2 = zeros(totalRows/rowsFromInd,1);

for ind = 1:size(rows,1)
    kaikkiRivit = rows(ind,1):rows(ind,2);
    for riviNumero = rowsFromInd:rowsFromInd:length(kaikkiRivit)
    %for riviNumero = rowsFromInd:rowsFromInd:length(rows{ind})
        %rivi = rows{ind}(riviNumero);
        rivi = kaikkiRivit(riviNumero);
        partitio2(rivi/rowsFromInd) = PARTITION(ind);
    end
end
PARTITION = partitio2;

%---------------------------------------------------------------


%--------------------------------------------------------------------


function newline = takeLine(description,width)
%Returns one line from the description: line ends to the first
%space after width:th mark.
newLine = description(1:width);
n = width+1;
while ~isspace(description(n)) & n<length(description)
    n = n+1;
end;
newline = description(1:n);


%-------------------------------------------------------------------


function changesInLogml=writeMixtureInfoPop(logml, rows, data, adjprior, ...
    priorTerm, outPutFile, inputFile, partitionSummary, popnames)

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global LOGDIFF;
ninds = size(rows,1);
npops =  size(COUNTS,3);
names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksilˆihin

if length(outPutFile)>0
    fid = fopen(outPutFile,'a');
else
    fid = -1;
    diary('baps4_output.baps'); % save in text anyway.
end

dispLine;
disp('RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:');
disp(['Data file: ' inputFile]);
disp(['Number of clustered groups: ' ownNum2Str(ninds)]);
disp(['Number of clusters in optimal partition: ' ownNum2Str(npops)]);
disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
disp(' ');
if (fid ~= -1)
    fprintf(fid,'%s \n', ['RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:']); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clustered groups: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Number of clusters in optimal partition: ' ownNum2Str(npops)]); fprintf(fid,'\n');
    fprintf(fid,'%s \n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
    fprintf(fid,'\n');
end

cluster_count = length(unique(PARTITION));
disp(['Best Partition: ']);
if (fid ~= -1)
    fprintf(fid,'%s \n',['Best Partition: ']); fprintf(fid,'\n');
end
for m=1:cluster_count
    indsInM = find(PARTITION==m);
    length_of_beginning = 11 + floor(log10(m));
    cluster_size = length(indsInM);
    
    if names
        text = ['Cluster ' num2str(m) ': {' char(popnames{indsInM(1)})];
        for k = 2:cluster_size
            text = [text ', ' char(popnames{indsInM(k)})];
        end;    
    else
        text = ['Cluster ' num2str(m) ': {' num2str(indsInM(1))];
        for k = 2:cluster_size
            text = [text ', ' num2str(indsInM(k))];
        end;
    end
    text = [text '}'];
    while length(text)>58
        %Take one line and display it.
        new_line = takeLine(text,58);
        text = text(length(new_line)+1:end);
        disp(new_line);
        if (fid ~= -1)
            fprintf(fid,'%s \n',[new_line]);
            fprintf(fid,'\n');
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
            fprintf(fid,'%s \n',[text]);
            fprintf(fid,'\n');
        end
    end;
end

if npops > 1
    disp(' ');
    disp(' ');
    disp('Changes in log(marginal likelihood) if group i is moved to cluster j:');
    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if group i is moved to cluster j:']); %fprintf(fid, '\n');
    end

    if names
        nameSizes = zeros(ninds,1);
        for i = 1:ninds
            nimi = char(popnames{i});
            nameSizes(i) = length(nimi);
        end
        maxSize = max(nameSizes);
        maxSize = max(maxSize, 5);
        erotus = maxSize - 5;
        alku = blanks(erotus);
        ekarivi = [alku 'group' blanks(6+erotus)];
    else
        ekarivi = 'group      ';
    end
    for i = 1:cluster_count
        ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
    end
    disp(ekarivi);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
    end

    changesInLogml = LOGDIFF';
    for ind = 1:ninds
        %[muutokset, diffInCounts] = laskeMuutokset(ind, rows, data, ...
        %    adjprior, priorTerm);
        muutokset = changesInLogml(:,ind);
        if names
            nimi = char(popnames{ind});
            rivi = [blanks(maxSize - length(nimi)) nimi ':'];
        else
            rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
        end
        for j = 1:npops
            rivi = [rivi '  ' logml2String(omaRound(muutokset(j)))];
        end
        disp(rivi);
        if (fid ~= -1)
            fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
        end
    end

    disp(' '); disp(' ');
    disp('KL-divergence matrix in PHYLIP format:');
    dist_mat = zeros(npops, npops);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['KL-divergence matrix in PHYLIP format:']); %fprintf(fid, '\n');
    end

    maxnoalle = size(COUNTS,1);
    nloci = size(COUNTS,2);
    d = zeros(maxnoalle, nloci, npops);
    prior = adjprior;
    prior(find(prior==1))=0;
    nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
    prior(1,nollia)=1;
    for pop1 = 1:npops
        d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
        %dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
    end
%     ekarivi = blanks(7);
%     for pop = 1:npops
%         ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
%     end
    ekarivi = num2str(npops);
    disp(ekarivi);

    if (fid ~= -1)
        fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
    end

    for pop1 = 1:npops
        rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
        for pop2 = 1:pop1-1
            dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
            div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
            div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
            div = (div12+div21)/2;
            % rivi = [rivi kldiv2str(div) '  '];
            dist_mat(pop1,pop2) = div;
        end
%         disp(rivi);
%         if (fid ~= -1)
%             fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
%         end
    end
   
else
    changesInLogml = [];
end

dist_mat = dist_mat + dist_mat'; % make it symmetric
for pop1 = 1:npops
    rivi = ['Cluster_' num2str(pop1) ' '];
    for pop2 = 1:npops
        rivi = [rivi kldiv2str(dist_mat(pop1,pop2)) ' '];
    end
    disp(rivi);
    if (fid ~= -1)
         fprintf(fid, '%s \n', [rivi]); %fprintf(fid, '\n');
    end
end

disp(' ');
disp(' ');
disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

if (fid ~= -1)
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
    fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values']); fprintf(fid, '\n');
end

partitionSummary = sortrows(partitionSummary,2);
partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
partitionSummary = partitionSummary(find(partitionSummary(:,2)>-1e49),:);
if size(partitionSummary,1)>10
    vikaPartitio = 10;
else
    vikaPartitio = size(partitionSummary,1);
end
for part = 1:vikaPartitio
    line = [num2str(partitionSummary(part,1)) '    ' num2str(partitionSummary(part,2))];
    disp(line);
    if (fid ~= -1)
        fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
    end
end

if ~fixedK

    disp(' ');
    disp(' ');
    disp('Probabilities for number of clusters');

    if (fid ~= -1)
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
        fprintf(fid, '%s \n', ['Probabilities for number of clusters']); fprintf(fid, '\n');
    end

    npopsTaulu = unique(partitionSummary(:,1));
    len = length(npopsTaulu);
    probs = zeros(len,1);
    partitionSummary(:,2) = partitionSummary(:,2)-max(partitionSummary(:,2));
    sumtn = sum(exp(partitionSummary(:,2)));
    for i=1:len
        npopstn = sum(exp(partitionSummary(find(partitionSummary(:,1)==npopsTaulu(i)),2)));
        probs(i) = npopstn / sumtn;
    end
    for i=1:len
        if probs(i)>1e-5
            line = [num2str(npopsTaulu(i)) '   ' num2str(probs(i))];
            disp(line);
            if (fid ~= -1)
                fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
            end
        end
    end
end

if (fid ~= -1)
    fclose(fid);
else
    diary off
end


%---------------------------------------------------------------


function dispLine;
disp('---------------------------------------------------');

%--------------------------------------------------------------

function num2 = omaRound(num)
% Pyˆrist‰‰ luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------

function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in?
% yks t‰ytyy olla kokonaisluku, joka on 
% v‰hint‰‰n -1:n suuruinen. Pienemmill?
% luvuilla tapahtuu jokin pyˆristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);


function mjono = kldiv2str(div)
mjono = '      ';
if abs(div)<100
    %Ei tarvita e-muotoa
    mjono(6) = num2str(rem(floor(div*1000),10));
    mjono(5) = num2str(rem(floor(div*100),10));
    mjono(4) = num2str(rem(floor(div*10),10));
    mjono(3) = '.';
    mjono(2) = num2str(rem(floor(div),10));
    arvo = rem(floor(div/10),10);
    if arvo>0
        mjono(1) = num2str(arvo);
    end
    
else
    suurinYks = floor(log10(div));
    mjono(6) = num2str(suurinYks);
    mjono(5) = 'e';
    mjono(4) = palautaYks(abs(div),suurinYks-1);
    mjono(3) = '.';
    mjono(2) = palautaYks(abs(div),suurinYks);
end


%-----------------------------------------------


function ninds = testaaOnkoKunnollinenBapsData(data)
%Tarkastaa onko viimeisess?sarakkeessa kaikki
%luvut 1,2,...,n johonkin n:‰‰n asti.
%Tarkastaa lis‰ksi, ett?on v‰hint‰‰n 2 saraketta.
if size(data,1)<2
    ninds = 0; return;
end
lastCol = data(:,end);
ninds = max(lastCol);
if ~isequal((1:ninds)',unique(lastCol))
    ninds = 0; return;
end

%--------------------------------------------------------------------------
%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

function initial_partition=admixture_initialization(data_matrix,nclusters, Z)

size_data=size(data_matrix);
nloci=size_data(2)-1;
n=max(data_matrix(:,end));
T=cluster_own(Z,nclusters);
initial_partition=zeros(size_data(1),1);
for i=1:n
    kori=T(i);
    here=find(data_matrix(:,end)==i);
    for j=1:length(here)
        initial_partition(here(j),1)=kori;
    end
end

function T = cluster_own(Z,nclust)
true=logical(1);
false=logical(0);
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

function [sumcounts, counts, logml] = ...
    initialCounts(partition, data, npops, rows, noalle, adjprior)

nloci=size(data,2);
ninds = size(rows, 1);

%koot = rows(:,1) - rows(:,2) + 1;
%maxSize = max(koot);

counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);
for i=1:npops
    for j=1:nloci
        havainnotLokuksessa = find(partition==i & data(:,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

%initializeGammaln(ninds, maxSize, max(noalle));

logml = laskeLoggis(counts, sumcounts, adjprior);

%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedet‰‰n, ett?annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss?ei viel?ole
% annettua logml arvoa, niin lis‰t‰‰n worstIndex:in kohtaan uusi logml ja
% nykyist?partitiota vastaava nclusters:in arvo. Muutoin ei tehd?mit‰‰n.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lˆydetty partitio ei ole viel?kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end

%--------------------------------------------------------------------------

function inds = returnInOrder(inds, pop, globalRows, data, ...
    adjprior, priorTerm)
% Palauttaa yksilˆt j‰rjestyksess?siten, ett?ensimm‰isen?on
% se, jonka poistaminen populaatiosta pop nostaisi logml:n
% arvoa eniten.

global COUNTS;      global SUMCOUNTS;
ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind =inds(i);
    rows = globalRows(i,1):globalRows(i,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop) = COUNTS(:,:,pop)-diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)-diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
    COUNTS(:,:,pop) = COUNTS(:,:,pop)+diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)+diffInSumCounts;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);

%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimm‰isen tyhj‰n populaation indeksin. Jos tyhji?
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;
pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end
