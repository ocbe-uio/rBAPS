function view_geneflow

%waitALittle;
% [filename, pathname] = uigetfile('*.mat', 'Load admixture results.');
% if (sum(filename)==0) || (sum(pathname)==0)
%     return;
% end

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
h0 = findobj('Tag','filename1_text');
filename = get(h0,'String');

% struct_array = load([pathname filename]);
disp('---------------------------------------------------');
disp('Viewing the gene flow.');
disp(['Load the admixture result from: ',[filename],'...']);
% if isfield(struct_array,'c')  %Matlab versio
%     c = struct_array.c;
%     if ~isfield(c,'proportionsIt')
%         disp('*** ERROR: Incorrect file format');
%         return
%     end
% elseif isfield(struct_array,'proportionsIt')  %Mideva versio
%     c = struct_array;
%     if ~isfield(c,'proportionsIt')
%         disp('*** ERROR: Incorrect file format');
%         return
%     end
% else
%     disp('*** ERROR: Incorrect file format');
%     return;
% end

proportionsIt = c.proportionsIt; 
popnames = c.popnames; partition = c.PARTITION;

% if isempty(popnames) || size(popnames,1)==size(partition,1) % bacteria/spatial    
    if isempty(popnames)
        ninds = size(partition,1);
        popnames=cell(ninds,2);
        for ind=1:ninds
            popnames{ind,1}=cellstr(num2str(ind));
        end        
        popnames(:,2)=num2cell((1:ninds)');
    end
    
    npops = c.npops;
    if isfield(c,'admixnpops') 
        admixnpops = c.admixnpops;
    else % if the admixture result is based on pre-defined partitions.
        admixnpops = npops;
    end
    
    if ~isfield(c,'pvalue') % compatiable with old data
        disp('*** WARNING: Old admixture format detected.');
        disp('*** WARNING: pvalue is not found in the admixture result.');
        disp('*** WARNING: all the admixture will be significant.');
        pvalue = ones(size(partition,1),1);
    else
        pvalue = c.pvalue;
    end
    
    removed_inds = logical(~any(proportionsIt,2));
    outliers = unique(partition(removed_inds));
    total_clusters = [1:npops];
    if ~isempty(outliers)
        fprintf('Outlier clusters: ');
        for i = 1:length(outliers)
            fprintf(['%s' blanks(2)], num2str(outliers(i)));
        end
        fprintf('\n');
        inliers = setdiff(total_clusters, outliers);
    else
        inliers = total_clusters;       
    end
    
    if length(inliers)~=admixnpops
        disp('*** ERROR: in view linkage admixture.');
        return
    end

    groupnames = cell(1,admixnpops);
    for i=1:admixnpops
        groupnames{i} = sprintf('Cluster %d',inliers(i));
    end

    %waitALittle;
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
                incluster = logical(proportionsIt(i,:)==max(proportionsIt(i,:)));
                proportionsIt(i,:) = zeros(1, admixnpops);
                proportionsIt(i,incluster) = 1;
            end
        end

        ix = find(ismember(partition(:),inliers));
        proportionsIt = proportionsIt(ix,:);
        partition = partition(ix);
    end
    
    % Calculate the gene flow weight.
    % W[j,i]: the gene flow from cluster j to i.
    W = ones(admixnpops,admixnpops); 
    for i=1:admixnpops
        inds = logical(partition==inliers(i));
        ninds = sum(inds);
        for j = 1:admixnpops
            W(j,i) = sum(proportionsIt(inds,j))/ninds;
        end
    end
    disp('Gene flow matrix (row=source col=target): ');
    ekarivi = 'Cluster      ';
    for i = 1:admixnpops
        ekarivi = [ekarivi ownNum2Str(inliers(i)) blanks(8-floor(log10(i)))];
    end
    disp(ekarivi);
    for i = 1:admixnpops
        rivi = [blanks(4-floor(log10(i))) num2str(inliers(i)) ':'];
        for j = 1:admixnpops
            rivi = [rivi '   ' num2str(omaRound(W(i,j)),'%5.4f')];
        end
        disp(rivi);
    end
    
    disp('Generating the gene flow graph. Please wait...');
    % Show graph 
%     %waitALittle;
%     talle = questdlg(['Do you want to view the gene flow graph?' ...
%         ], 'Visualization?', 'Yes', 'No', 'Yes');
%     if isequal(talle,'No')
%         return
%     else
        h0 = findobj('Tag','geneflow_menu');
        graphviz_path = get(h0,'Userdata');
        if isempty(graphviz_path)
            %waitALittle;
            graphviz_path = uigetdir('','Specify the location of dot.exe in the GraphViz package: ');
            if graphviz_path == 0
                return;
            else
                h0 = findobj('Tag','geneflow_menu');
                set(h0,'Userdata',graphviz_path);
            end
        end
            % store the current directory.
            d = cd;
            try
                graphvis2(W,filename,inliers,graphviz_path,npops);
                disp('Finished.');
            catch
                disp('*** ERROR: dot.exe was not found.');
                h_errdlg = errordlg('dot.exe was not found in the specified path.');
                set(h_errdlg,'WindowStyle','Modal')
                set(h0,'Userdata',[]);
            end
            cd(d);
%     end
% else 
%    disp('This module is under developpment.');
% end

% -------------------------------------------------------------------------
function num2 = omaRound(num)
% Pyöristää luvun num 1 desimaalin tarkkuuteen
num = num*10000;
num = round(num);
num2 = num/10000;




