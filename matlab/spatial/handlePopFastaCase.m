function handlePopFastaCase(cc,pgPart,pgDist)
% specicially written to handle FASTA file format of individual clustering
% Lu Cheng, 15.12.2012

OUTPUT_FILE = 'baps6_output.baps';

teksti = 'Input upper bound to the number of populations (only one value): ';
npopstextExtra = inputdlg(teksti ,'Input maximum number of populations',1,{'20'});
if isempty(npopstextExtra)  % Painettu Cancel:ia
    return
else
    nMaxPops = str2num(npopstextExtra{1});
    nMaxPops = nMaxPops(1);
end

nPregroup = length(unique(pgPart));

roundTypes = [2*ones(1,nMaxPops) ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 1 1 1 1 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];

[partition, logml, partitionSummary, logmldiff] = model_search_pregroup(cc, pgPart, pgDist, roundTypes, nMaxPops);

cc.PARTITION = partition;  %note that the partition only contain nPregroup elements
cc.npops = length(unique(partition));
cc.logml = logml;
cc.partitionSummary = partitionSummary;
cc.logmldiff = logmldiff;

if cc.npops==nMaxPops
    choice = questdlg(sprintf('%d populations discovered, which is the same as input. We suggest you to set a larger number. Do you want to quit?', cc.npops),...
        'Yes''No','Yes');
    if strcmp(choice,'Yes')
        return
    end
end

writeMixtureInfo(cc);

popnames = cc.popnames;
pointers = cc.pointers;
vorPoints = cc.vorPoints;
vorCells = cc.vorCells;
coordinates = cc.coordinates;

if isequal(popnames, [])
    names = pointers;
else
   names = cell(size(pointers));
   indices = zeros(size(popnames(:,2)));
   for i=1:length(popnames(:,2));
       indices(i) = popnames{i,2};
   end
   for i = 1:length(pointers)
       inds = pointers{i};
       namesInCell = [];
       for j = 1:length(inds)
           ind = inds(j);
           I = find(indices > ind);
           if isempty(I)
               nameIndex = length(indices);
           else
               nameIndex = min(I) -1;
           end
           name = popnames{nameIndex};
           namesInCell = [namesInCell name];
       end
       names{i} = namesInCell;
   end
end
viewMixPartition(partition, popnames);
vorPlot(vorPoints, vorCells, partition, pointers, coordinates, names);    

talle = questdlg(['Do you want to save the mixture populations ' ...
    'so that you can use them later in admixture analysis or plot ' ...
    'additional images?'], ...
    'Save results?','Yes','No','Yes');
if isequal(talle,'Yes')
    %%waitALittle;    % Hetki odotusta, jotta muistaa kysy?..
    [filename, pathname] = uiputfile('*.mat','Save results as');
    
    if (filename == 0) & (pathname == 0)
        % Cancel was pressed
        return
    else % copy 'baps4_output.baps' into the text file with the same name.
        if exist(OUTPUT_FILE,'file')
            copyfile(OUTPUT_FILE,[pathname filename '.txt'])
            delete(OUTPUT_FILE)
        end
    end
       
    %  added by Lu Cheng, 05.12.2012
    tmpFile = [pathname filename '.mapfile.txt'];
    fid = fopen(tmpFile,'w+');
    fprintf(fid,'GroupLabel\tLatitude\tLongitude\tDescription\tLabel\n');
    for i=1:nPregroup
        fprintf(fid,'%d\t%.10f\t%.10f\t%d_%d\t%d\n',i,coordinates(i,1),coordinates(i,2),...
                i,partition(i),partition(i));
    end
    fclose(fid);
    
%     save([pathname filename], 'c');
    format_type = 'FASTA';
    save([pathname filename], 'cc','partition','pgDist','pgPart','format_type','-v7.3');
else
    if exist(OUTPUT_FILE,'file')
        delete(OUTPUT_FILE)
    end
end



%%%%%%%%%%%%%
function writeMixtureInfo(c)

outputFile = 'baps6_output.baps';

% output the semi-supervised clustering results to the outputFile
% modified by Lu Cheng, 28.03.2010

ninds = length(c.PARTITION);
npops =  c.npops;
popnames = c.popnames;
logml = c.logml;
partition = c.PARTITION;
partitionSummary = c.partitionSummary;

if isempty(popnames)
    popnames = cell(c.nPregroup,1);
    for i=1:c.nPregroup
        popnames{i} = num2str(i);
    end
end

if ~isempty(outputFile)
    fid = fopen(outputFile,'w+');
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
    indsInM = unique(c.groupPartition(partition==m));
    
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

logmldiff = c.logmldiff;
if npops == 1
    logmldiff = [];
else
    disp(' ');
    disp(' ');
    disp('Changes in log(marginal likelihood) if pregroup i is moved to cluster j:');
        
    if (fid ~= -1)
        fprintf(fid, '%s \n', ' '); fprintf(fid, '\n');
        fprintf(fid, '%s \n', 'Changes in log(marginal likelihood) if indvidual i is moved to cluster j:'); fprintf(fid, '\n');
    end

    text = sprintf('%10s','ind');
    for ii = 1:npops
        tmpstr = sprintf('\t%10s',num2str(ii));
        text = [text tmpstr];
    end
    
    disp(text);
    if (fid ~= -1)
        fprintf(fid, '%s \n', text);
    end
        
    for ii = 1:c.nPregroup
        text = sprintf('%10s',popnames{ii});
        for jj = 1:npops
            tmpstr = sprintf('\t%10s',num2str(logmldiff(ii,jj),'%10.6f'));
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

partitionSummaryKaikki = partitionSummary;
partitionSummary =[];
for i=1:size(partitionSummaryKaikki,3)
    partitionSummary = [partitionSummary; partitionSummaryKaikki(:,:,i)];
end
% [I,J] = find(partitionSummaryKaikki(:,2,:)>-1e49);
% partitionSummaryKaikki = partitionSummaryKaikki(I,:,:);

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


%%%%%%%%%%


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

