function admix2

global PARTITION; global COUNTS;
global SUMCOUNTS;
clearGlobalVars;

input_type = questdlg('Specify the format of your data: ',...
    'Specify Data Format', ...
    'BAPS-format', 'GenePop-format', 'BAPS-format');

switch input_type
   
case 'BAPS-format'
    waitALittle;
    [filename, pathname] = uigetfile('*.txt', 'Load data in BAPS-format');
    if filename==0
        return;
    end
	
    data = load([pathname filename]);
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    if (ninds==0) 
        disp('Incorrect Data-file.');
        return;    
    end
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;
    waitALittle;
    [filename, pathname] = uigetfile('*.txt', 'Load Partition');
    if filename==0
        return;
    end
    PARTITION = load([pathname filename]);
    if ~(size(PARTITION,2)==1) | ~(size(PARTITION,1)==ninds)
        disp('Incorrect Partition-file.');
        return;
    end
    
    input_pops = questdlg(['When using data which are in BAPS-format, '...
       'you can specify the sampling populations of the individuals by '...
       'giving two additional files: one containing the names of the '...
       'populations, the other containing the indices of the first '...
       'individuals of the populations. Do you wish to specify the '...
       'sampling populations?'], ...
       'Specify sampling populations?',...
       'Yes', 'No', 'No');
    if isequal(input_pops,'Yes')
        waitALittle;
        [namefile, namepath] = uigetfile('*.txt', 'Load population names');
        if namefile==0
            kysyToinen = 0;
        else
            kysyToinen = 1;
        end
        if kysyToinen==1
            waitALittle;
            [indicesfile, indicespath] = uigetfile('*.txt', 'Load population indices');
            if indicesfile==0
                popnames = [];     
            else
                popnames = initPopNames([namepath namefile],[indicespath indicesfile]);
            end
        else
            popnames = [];
        end                
    else
       	popnames = [];
    end
	   
    [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
    data = data(:, 1:end-1);
    npops = length(unique(PARTITION(find(PARTITION>=0))));
    
    
case 'GenePop-format'
    waitALittle;
    [filename, pathname] = uigetfile('*.txt', 'Load data in GenePop-format');
    if filename==0
        return;
    end
    
    kunnossa = testaaGenePopData([pathname filename]);
    if kunnossa==0
        return
    end

    [data,popnames]=lueGenePopData([pathname filename]);
	
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;
    
    [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
    data = data(:, 1:end-1);
    
    npops = size(popnames,1);
    ninds = size(data,1)/rowsFromInd;
    PARTITION = zeros(ninds,1);
    ind = 1;
    for pop = 1:npops-1
        while (ind < popnames{pop+1,2})
            PARTITION(ind) = pop;
            ind = ind+1;
        end
    end
    while (ind <= ninds)
        PARTITION(ind) = npops;
        ind = ind+1;
    end

    all_in_text = questdlg(['Do you wish to use also the last population in the ',...
        'data to define one more population for admixture analysis: '],...
        'Define a population based on the last population in the data?', ...
        'Yes', 'No', 'Yes');
    if isequal(all_in_text, 'No')
        PARTITION(find(PARTITION==npops)) = -1;
        npops = npops-1;
    end
    otherwise return
end

initialPartition = PARTITION(:,ones(1,rowsFromInd))';
initialPartition = initialPartition(:);
[sumcounts, counts, logml] = ...
    initialCounts(initialPartition, data, npops, rowsFromInd, noalle);
COUNTS = counts; SUMCOUNTS = sumcounts;

clear('initialPartition', 'counts', 'sumcounts', ... 
    'filename', 'ind', 'input_type', ... 
    'logml', 'ninds', 'pathname', 'pop', 'priorTerm');
clear('indicesfile','indicespath','input_pops','kysyToinen',...
    'namefile','namepath');
c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
c.alleleCodes = alleleCodes; c.adjprior = adjprior; c.popnames = popnames;
c.rowsFromInd = rowsFromInd; c.data = data; c.npops = npops; c.noalle = noalle;
admix1(c);

%------------------------------------------------------------------------


function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];

%--------------------------------------------------------

function ninds = testaaOnkoKunnollinenBapsData(data)
%Tarkastaa onko viimeisessä sarakkeessa kaikki
%luvut 1,2,...,n johonkin n:ään asti.
%Tarkastaa lisäksi, että on vähintään 2 saraketta.
if size(data,1)<2
    ninds = 0; return;
end
lastCol = data(:,end);
ninds = max(lastCol);
if ~isequal((1:ninds)',unique(lastCol))
    ninds = 0; return;
end

%-----------------------------------------------------------------------------------


function popnames = initPopNames(nameFile, indexFile)
%Palauttaa tyhjän, mikäli nimitiedosto ja indeksitiedosto
% eivät olleet yhtä pitkiä.

popnames = [];
indices = load(indexFile);

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

if length(names) ~= length(indices)
    disp('The number of population names must be equal to the number of ');
    disp('entries in the file specifying indices of the first individuals of ');
    disp('each population.');
    return;
end

popnames = cell(length(names), 2);
for i = 1:length(names)
    popnames{i,1} = names(i);
    popnames{i,2} = indices(i);
end

%---------------------------------------------------------------------------------------


function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
    handleData(raw_data)
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

adjprior = zeros(max(noalle),nloci);
priorTerm = 0;
for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j)); 
end

%--------------------------------------------------------------------


function kunnossa = testaaGenePopData(tiedostonNimi)
% kunnossa == 0, jos data ei ole kelvollinen genePop data.
% Muussa tapauksessa kunnossa == 1.

kunnossa = 0;
fid = fopen(tiedostonNimi);
line1 = fgetl(fid);  %ensimmäinen rivi
line2 = fgetl(fid);  %toinen rivi
line3 = fgetl(fid);  %kolmas

if (isequal(line1,-1) | isequal(line2,-1) | isequal(line3,-1))
    disp('Incorrect file format 1168'); fclose(fid);
    return 
end
if (testaaPop(line1)==1 | testaaPop(line2)==1)
    disp('Incorrect file format 1172'); fclose(fid);
    return
end
if testaaPop(line3)==1
    %2 rivi tällöin lokusrivi
    nloci = rivinSisaltamienMjonojenLkm(line2);
    line4 = fgetl(fid);
    if isequal(line4,-1)
        disp('Incorrect file format 1180'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin neljä täytyy sisältää pilkku.
        disp('Incorrect file format 1185'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetään, että pysähtyy
        pointer = pointer+1;
    end
    line4 = line4(pointer+1:end);  %pilkun jälkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format 1195'); fclose(fid);
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
        disp('Incorrect file format 1206'); fclose(fid);
        return
    end
    nloci = lineNumb-2;
    line4 = fgetl(fid);  %Eka rivi pop sanan jälkeen
    if isequal(line4,-1)
        disp('Incorrect file format 1212'); fclose(fid);
        return 
    end
    if ~any(line4==',')
        % Rivin täytyy sisältää pilkku.
        disp('Incorrect file format 1217'); fclose(fid);
        return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetään, että pysähtyy.
        pointer = pointer+1;
    end
 
    line4 = line4(pointer+1:end);  %pilkun jälkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
        disp('Incorrect file format 1228'); fclose(fid);
        return
    end
end
kunnossa = 1;
fclose(fid);

%------------------------------------------------------
    

function [data, popnames] = lueGenePopData(tiedostonNimi)

fid = fopen(tiedostonNimi);
line = fgetl(fid);  %ensimmäinen rivi
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
        popnames{nimienLkm, 1} = {nimi};   %Näin se on greedyMix:issäkin?!?
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

data = data(1:ninds*2,:);
popnames = popnames(1:nimienLkm,:);
fclose(fid);


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
% line on ensimmäinen pop-sanan jälkeinen rivi
% Genepop-formaatissa olevasta datasta. funktio selvittää
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
% Palauttaa line:n sisältämien mjonojen lukumäärän.
% Mjonojen välissä täytyy olla välilyönti.
count = 0;
pit = length(line);
tila = 0;    %0, jos odotetaan välilyöntejä, 1 jos odotetaan muita merkkejä
for i=1:pit
    merkki = line(i);
    if (isspace(merkki) & tila==0) 
        %Ei tehdä mitään.
    elseif (isspace(merkki) & tila==1)
        tila = 0;
    elseif (~isspace(merkki) & tila==0)
        tila = 1;
        count = count+1;
    elseif (~isspace(merkki) & tila==1)
        %Ei tehdä mitään
    end
end

%-------------------------------------------------------

function pal = testaaPop(rivi)
% pal=1, mikäli rivi alkaa jollain seuraavista
% kirjainyhdistelmistä: Pop, pop, POP. Kaikissa muissa
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

%-----------------------------------------------------------------------------------


function [sumcounts, counts, logml] = ...
    initialCounts(partition, data, npops, rowsFromInd, noalle)

nloci=size(data,2);
ninds = size(data,1)/rowsFromInd;

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

initializeGammaln(ninds, rowsFromInd, max(noalle));

logml = computeLogml(counts, sumcounts, noalle, data, rowsFromInd);


%-----------------------------------------------------------------------


function logml=computeLogml(counts, sumcounts, noalle, data, rowsFromInd)
nloci = size(counts,2);
npops = size(counts,3);
adjnoalle = zeros(max(noalle),nloci);
for j=1:nloci
    adjnoalle(1:noalle(j),j)=noalle(j);
    if (noalle(j)<max(noalle))
        adjnoalle(noalle(j)+1:end,j)=1;
	end       
end

%logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
%    - npops*sum(sum(gammaln(adjprior))) - ...
%    sum(sum(gammaln(1+sumcounts)));
%logml = logml2;

global GAMMA_LN;
rowsInG = size(data,1)+rowsFromInd;

logml = sum(sum(sum(GAMMA_LN(counts+1 + repmat(rowsInG*(adjnoalle-1),[1 1 npops]))))) ...
    - npops*sum(sum(GAMMA_LN(1, adjnoalle))) ...
    -sum(sum(GAMMA_LN(sumcounts+1,1)));


%--------------------------------------------------------------------------


function initializeGammaln(ninds, rowsFromInd, maxAlleles)
%Alustaa GAMMALN muuttujan s.e. GAMMALN(i,j)=gammaln((i-1) + 1/j)
global GAMMA_LN;
GAMMA_LN = zeros((1+ninds)*rowsFromInd, maxAlleles);
for i=1:(ninds+1)*rowsFromInd
    for j=1:maxAlleles
        GAMMA_LN(i,j)=gammaln((i-1) + 1/j);
    end
end

%--------------------------------------------------------


function data = addAlleles(data, ind, line, divider)
% Lisaa BAPS-formaatissa olevaan datataulukkoon
% yksilöä ind vastaavat rivit. Yksilön alleelit
% luetaan genepop-formaatissa olevasta rivistä
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
