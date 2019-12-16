function view_admixture(proportionsIt,npops,admixnpops,...
    popnames,partition,pvalue,source)

% Find the outlier individuals
% removed_inds = find(~any(proportionsIt,2));
removed_inds = logical(~any(proportionsIt,2));
outliers = unique(partition(removed_inds));
if ~isempty(outliers)
    fprintf('Outlier clusters: ');
    for i = 1:length(outliers)
        fprintf(['%s' blanks(2)], num2str(outliers(i)));
    end
    fprintf('\n');
end

% total_clusters = [1:max(partition)];
total_clusters = [1:npops];
inliers = setdiff(total_clusters, outliers);
if length(inliers)~=admixnpops
    disp('*** ERROR: in view linkage admixture.');
    return
end

if isempty(popnames) || size(popnames,1)==size(partition,1)
    
    groupnames = cell(1,admixnpops);
    for i=1:admixnpops
        groupnames{i} = sprintf('Cluster %d',inliers(i));
    end
    waitALittle;
    [s,v] = listdlg('PromptString',[sprintf('%d ',admixnpops) 'Clusters available:'],...
        'SelectionMode','multiple',...
        'Name','Select Clusters',...
        'ListString',groupnames);
    if isempty(s) || ~v
        disp('*** WARNING: Viewing admixture cancelled.');
        return
    else
        fprintf('Viewing admixture in cluster(s): %s\n',num2str(inliers(s)));
    end

    waitALittle;
    answer = inputdlg( ['Input the p value for the strain admixture. Strains with '...
        'nonsignificant admixture will be displayed as a single color bar'],...
        'Input the upper bound p value',1,{'0.05'});
    if isempty(answer)  % cancel has been pressed
        return
    else
        maxp = str2num(answer{1});
        fprintf('P value: %4.2f\n', maxp);

        % nonsignificant strains are forced to be single color bar.
        nonsignificant_ix = find( (pvalue>maxp & pvalue <=1) );
        if ~isempty(nonsignificant_ix)
            for i = nonsignificant_ix'
                % incluster = find(proportionsIt(i,:)==max(proportionsIt(i,:)));
                incluster = logical(proportionsIt(i,:)==max(proportionsIt(i,:)) & proportionsIt(i,:)~=0);
                proportionsIt(i,:) = zeros(1, admixnpops);
                proportionsIt(i,incluster) = 1;
            end
        end

        ix = find(ismember(partition(:),inliers(s)));
        %         if isempty(ix)
        %             disp('*** ERROR: No significant data. Try again with a higher p value.');
        %             return
        %         end
        proportionsIt_raw = proportionsIt; 
        partition_raw = partition;
        proportionsIt = proportionsIt(ix,:);
        if isempty(popnames) || size(popnames,1)==size(partition,1)
            popnames = popnames(ix,:);
        end
        partition = partition(ix);
        pvalue = pvalue(ix);
    end
else
    waitALittle;
    answer = inputdlg( ['Input the p value for the strain admixture. Strains with '...
        'nonsignificant admixture will be displayed as a single color bar'],...
        'Input the upper bound p value',1,{'0.05'});
    if isempty(answer)  % cancel has been pressed
        return
    else
        maxp = str2num(answer{1});
        fprintf('P value: %4.2f\n', maxp);

        % nonsignificant strains are forced to be single color bar.
        nonsignificant_ix = find( (pvalue>maxp & pvalue <=1) );
        if ~isempty(nonsignificant_ix)
            for i = nonsignificant_ix'
                % incluster = find(proportionsIt(i,:)==max(proportionsIt(i,:)));
                incluster = logical(proportionsIt(i,:)==max(proportionsIt(i,:)) & proportionsIt(i,:)~=0);
                proportionsIt(i,:) = zeros(1, admixnpops);
                proportionsIt(i,incluster) = 1;
            end
        end

        % ix = find(ismember(partition(:),inliers(s)));
        %         if isempty(ix)
        %             disp('*** ERROR: No significant data. Try again with a higher p value.');
        %             return
        %         end
        % proportionsIt_raw = proportionsIt; 
        % partition_raw = partition;
        % proportionsIt = proportionsIt(ix,:);
%         if isempty(popnames) || size(popnames,1)==size(partition,1)
%             popnames = popnames(ix,:);
%         end
%         partition = partition(ix);
%         pvalue = pvalue(ix);
    end

end



talle = questdlg(['Do you want names to be visible in the admixture ' ...
    'result graphics?'], 'Names visible?', 'Yes', 'No', 'Yes');

if isequal(talle,'No')
    if isempty(popnames) || size(popnames,1)==size(partition,1)
        viewPartition4(proportionsIt, [], npops, ...
            admixnpops, inliers, partition, source);
    else
        viewPartition2(proportionsIt, [], admixnpops, partition, source);
    end
else
    if isempty(popnames) || size(popnames,1)==size(partition,1)
        viewPartition4(proportionsIt, popnames, npops, ...
            admixnpops, inliers, partition, source);
    else
        viewPartition2(proportionsIt, popnames, admixnpops, partition, source);
    end
end

if isempty(popnames) || size(popnames,1)==size(partition,1)
    talle = questdlg(['Permutate the labelling orders to get a more structured graphics? '], 'Reorder the sample?', 'Yes', 'No', 'Yes');

    if isequal(talle,'Yes')
        h0 = findobj('Tag','admixture_result');
        if ~isempty(h0)
            close(h0);
            % drawnow;
            waitALittle;
            [partition_ordered, proportionsIt,popnames,pvalue] = ...
                viewPartition3(proportionsIt, popnames, npops, admixnpops, ...
                inliers, partition, pvalue, source);
            %         c.mixtureType = 'linkage_mix';
            %         c.npops = npops;
            talle = questdlg('Save the ordered admixture result in text?', ...
                'Save results?','Yes','No','Yes');
            if isequal(talle,'Yes')
                waitALittle;
                [filename, pathname] = uiputfile('*.txt','Save results as');
                if isempty(filename) && isempty(pathname)
                    % Cancel was pressed
                    return
                else
                    fprintf(1,'Saving the admixture result...\n');
                    dest = [pathname filename];
                    clusters = inliers(s);
                    tulostaAdmixtureTiedot(popnames, proportionsIt, pvalue, source, dest,...
                        maxp,clusters);
                end
            end

        else
            disp('***ERROR: In plotting the admixture results.');
            return;
        end
        drawnow;
    end
end


%--------------------------------------------------------------------------


function tulostaAdmixtureTiedot(popnames, proportions, uskottavuus, ...
                                  source, dest, maxp,clusters)

fid = fopen(dest,'w');
ninds = length(uskottavuus);
npops = size(proportions,2);

if fid ~= -1
    fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['--------------------------------------------']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['RESULTS OF ADMIXTURE ANALYSIS BASED']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['ON MIXTURE CLUSTERING OF INDIVIDUALS']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Source file: ' source]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Total number of clusters: ' npops]); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Viewing admixture in cluster(s): ',num2str(clusters)]);fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Number of individuals in selection: ' num2str(ninds)]); fprintf(fid, '\n');
    fprintf(fid,'P value: %4.2f \n', maxp); fprintf(fid, '\n');
    fprintf(fid, '\n');
end

namelength = zeros(1,ninds);
for ind = 1:ninds
    namelength(ind) = length(popnames{ind,1}{1});
end
maxlength = max(namelength);
ekaRivi = ['Index' blanks(2) blanks(maxlength-4) 'Name' blanks(2)];
for pop = 1:npops
    ekaRivi = [ekaRivi blanks(3-floor(log10(pop))) num2str(pop) blanks(2)];
end
ekaRivi = [ekaRivi blanks(2) 'p']; % Added on 29.08.06
disp(ekaRivi);
fprintf(fid, '%s \n',ekaRivi); fprintf(fid,'\n');
for ind = 1:ninds
    index = num2str(popnames{ind,2});
    namelength=['%' num2str(2+maxlength) 's'];
    rivi = [sprintf('%5s',index) sprintf(namelength,popnames{ind,1}{1})...
                       blanks(2)];
    if any(proportions(ind,:)>0)
        for pop = 1:npops-1
            rivi = [rivi proportion2str(proportions(ind,pop)) blanks(2)];
        end
        rivi = [rivi proportion2str(proportions(ind,npops))];
        rivi = [rivi blanks(2) ownNum2Str(uskottavuus(ind))];
    end
    disp(rivi);
    if fid ~= -1
        fprintf(fid,'%s \n',rivi); fprintf(fid,'\n');
    end
end
fclose(fid);
fprintf(1,'finished.\n');

%--------------------------------------------------------------------------
% function dispLine
% disp('---------------------------------------------------');

%--------------------------------------------------------------------------
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


