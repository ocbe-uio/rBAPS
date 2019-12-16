function view_energy()

%waitALittle;
% [filename1, pathname1] = uigetfile('*.mat', 'Load mixture results.');
% if (sum(filename1)==0) || (sum(pathname1)==0)
%     return;
% end
% struct_array = load([pathname1 filename1]);

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
h0 = findobj('Tag','filename1_text');
filename1 = get(h0,'String');
disp('---------------------------------------------------');
disp('Viewing the energy lanscape (Kolmogorov–Smirnov test).');
disp(['Load the mixture result from: ',[filename1],'...']);

 if isfield(c,'c')  %Matlab versio
     c = c.c;
     if ~isfield(c,'mixtureType')
         disp('*** ERROR: Incorrect file format');
         return
     end
 %else
 %    disp('*** ERROR: Incorrect file format');
 %    return;
 end
if ~isfield(c,'changesInLogml')
    disp('*** WARNING: Old mixture format detected.');
    if 0  % Disabled 12/06/07, Jukka Siren
        switch c.mixtureType
            case 'linkage_mix'
                %waitALittle;
                [filename, pathname] = uigetfile('*.mat', 'Load the corresponding preprocessed data.');
                if (sum(filename)==0) || (sum(pathname)==0)
                    return;
                end
                struct_array = load([pathname filename]);
                disp('The corresponding preprocessed data is needed.');
                disp(['Load the preprocessed data from: ',[pathname filename], '...']);
                if isfield(struct_array,'c')
                    c_raw = struct_array.c;
                    if ~all(size(c_raw.adjprior)==size(c.adjprior)) ||...
                            ~all(size(c_raw.data)==size(c.data))
                        disp('Mixture result and preprocessed data do not match.');
                        return;
                    end

                    if ~isfield(c_raw,'linkage_model')
                        model_entry = input('Specify the linkage model(1-linear 2-codon):');
                        switch model_entry,
                            case 1
                                linkage_model = 'linear';
                            case 2
                                linkage_model = 'codon';
                        end
                    else
                        linkage_model = c_raw.linkage_model;
                    end
                else
                    disp('*** ERROR: Incorrect file format');
                    return;
                end
                disp('Start calculating log likelihood. Please wait...');
                npops = c.npops;
                partition = c.PARTITION;
                cq_counts = c.CQ_COUNTS;
                cq_sumcounts = c.CQ_SUMCOUNTS;
                sp_counts = c.SP_COUNTS;
                sp_sumcounts = c.SP_SUMCOUNTS;

                % Shrink the raw data
                noalle_cq = size(c.adjprior_cq,1);
                informative_cq = logical(sum(c.adjprior_cq)~=noalle_cq);
                noalle_sp = size(c.adjprior_sp,1);
                informative_sp = logical(sum(c.adjprior_sp)~=noalle_sp);
                cq_counts = uint16(cq_counts(:,informative_cq,:));
                cq_sumcounts = uint16(cq_sumcounts(:,informative_cq));
                sp_counts = uint16(sp_counts(:,informative_sp,:));
                sp_sumcounts = uint16(sp_sumcounts(:,informative_sp));

                data = c_raw.data;
                ninds = size(data,1);
                component_mat = c_raw.component_mat;
                clear c_raw;
                index = data(:,end);
                [data_clique, data_separator, noalle_clique, noalle_separator] = ...
                    transform4(data, component_mat, linkage_model);

                if length(noalle_clique)~=length(find(informative_cq))
                    disp('*** ERROR: The linkage model is not consistent with the data.');
                    return
                end

                data_clique = [data_clique index];
                data_separator = [data_separator index];

                [counts_cq, nalleles_cq, prior_cq, adjprior_cq]...
                    = allfreqsnew2(data_clique, double(noalle_clique));
                clear data_clique;
                [counts_sp, nalleles_sp, prior_sp, adjprior_sp]...
                    = allfreqsnew2(data_separator, double(noalle_separator));
                clear data_separator;
                counts_cq = uint8(counts_cq);
                counts_sp = uint8(counts_sp);

                pop_logml = computePopulationLogml(double(cq_counts), double(cq_sumcounts),...
                    double(sp_counts), double(sp_sumcounts),...
                    [1:npops], adjprior_cq, adjprior_sp);

                changesInLogml = zeros(npops,ninds);
                for ind = 1:ninds
                    indCqCounts = uint16(counts_cq(:,:,ind));
                    indSpCounts = uint16(counts_sp(:,:,ind));
                    changesInLogml(:,ind) = computeChanges(cq_counts, cq_sumcounts,...
                        sp_counts, sp_sumcounts,...
                        partition, pop_logml,...
                        ind, adjprior_cq, ...
                        adjprior_sp, indCqCounts, indSpCounts);
                end
            case 'mix' % Independence model.
                %waitALittle;
                [filename, pathname] = uigetfile('*.mat', 'Load the corresponding preprocessed data.');
                if (sum(filename)==0) || (sum(pathname)==0)
                    return;
                end
                struct_array = load([pathname filename]);
                disp('The corresponding preprocessed data is needed.');
                disp(['Load the preprocessed data from: ',[pathname filename], '...']);
                if isfield(struct_array,'c')
                    c_raw = struct_array.c;
                    if ~all(size(c_raw.adjprior)==size(c.adjprior)) ||...
                            ~all(size(c_raw.data(:,[1:end-1]))==size(c.data))
                        disp('Mixture result and preprocessed data do not match.');
                        return;
                    end
                else
                    disp('*** ERROR: Incorrect file format');
                    return;
                end


                npops = c.npops;
                pops = 1:npops;
                rowsFromInd = c.rowsFromInd;
                data = c.data;
                ninds = size(data,1)/rowsFromInd;
                adjprior = c.adjprior;
                priorTerm = c_raw.priorTerm;
                COUNTS = c.COUNTS;
                SUMCOUNTS = c.SUMCOUNTS;
                PARTITION = c.PARTITION;
                POP_LOGML = computePopulationLogml2(pops, adjprior, priorTerm, ...
                    COUNTS, SUMCOUNTS);
                changesInLogml = zeros(npops,ninds);
                for ind = 1:ninds
                    [muutokset, diffInCounts] = laskeMuutokset(ind, rowsFromInd, data, ...
                        adjprior, priorTerm, COUNTS, SUMCOUNTS, PARTITION, POP_LOGML);
                    changesInLogml(:,ind) = muutokset;
                end
            otherwise
                disp('This model is under construction.');
                return
        end


        c.changesInLogml = changesInLogml;

        fprintf(1,'Saving the result...')
        %waitALittle;
        save_preproc = questdlg('Do you wish to save the updated mixture result?',...
            'Save mixture results?',...
            'Yes','No','Yes');
        if isequal(save_preproc,'Yes');
            %waitALittle;
            [filename, pathname] = uiputfile('*.mat','Save mixture result as');
            if isempty(filename) && isempty(pathname)
                return;
            else
                kokonimi = [pathname filename];
%                 save(kokonimi,'c');
                save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
                filename1 = filename;
            end

        end;
        fprintf(1,'Finished.\n');
    end
    return
end

if isequal(c.mixtureType,'spatialPop') || isequal(c.mixtureType,'popMix')
    % Group level clustering, using the correct partition
    view_density(c.changesInLogml, c.groupPartition, filename1);
else
    view_density(c.changesInLogml, c.PARTITION, filename1);
end


% -------------------------------------------------------------------------
function view_density(changesInLogml, partition, filename)
npops = size(changesInLogml,1);
groupnames = cell(1,npops+1);
for i=1:npops
    groupnames{i} = sprintf('Cluster %d',i);
end

disp('Genetic affinity matrix (row=source col=target): test statistic(p value) ');
ekarivi = 'Cluster      ';
for i = 1:npops
    ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
end
disp(ekarivi);

    
for i=1:npops
    indsi = logical(partition==i); % choose cluster i
    nindsi = sum(indsi);
    rivi = [blanks(4-floor(log10(i))) num2str(i) ':'];
    
    h = zeros(npops,1); % indicator 0 same 1 different
    h1 = h;
    p = zeros(npops,1); % p values
    k = zeros(npops,1); % test statistic
    k1 = k;
    p2 = zeros(npops,1);
    for j=1:npops
        indsj = logical(partition==j);
        [h(j), p(j), k(j)] = kstest2(zscore(changesInLogml(j,indsi)), zscore(changesInLogml(i,indsj)));
        [h1(j), p1(j), k1(j)] = kstest2(zscore(changesInLogml(j,indsi)), zscore(changesInLogml(i,indsj)),0.05,'larger');
        [h2(j), p2(j), k2(j)] = kstest2(zscore(changesInLogml(j,indsi)), zscore(changesInLogml(i,indsj)),0.05,'smaller');


        if k(j)==k2(j)
            k(j) = -k(j); 
       
        end
        % [h2(j),p2(j),k2(j)] = kstest2(changesInLogml(j,indsi), changesInLogml(i,indsj),0.05,'smaller');
        % nonzero(i,j) = k1(j)+k2(j);
        % nonzero(i,j) =
        % (k(j)*(h1(j)-0.5))*(mean(changesInLogml(j,indsi)-mean(changesInLogml(i,indsj))));
        nonzero(i,j) = k(j);
        rivi = [rivi '   ' num2str(nonzero(i,j),'%5.4f') '(' num2str(p(j),'%5.4f') ')'];
    end
    disp(rivi);
    
    % sibling = setdiff(find(h==0),i);
    
    % parent = find(p==min(p));
%     fprintf(1,'Cluster %d           ',i);
%     fprintf(1,'%s(%f)            ',mat2str(parent),min(p));
%     
%     offspring = find(p2==min(p2));
%     fprintf(1,'%s(%f)\n',mat2str(offspring), min(p2));
end

% groupnames{npops+1} = sprintf('All clusters');
% %waitALittle;
% [s1,v1] = listdlg('PromptString','Select one source cluster:',...
%     'SelectionMode','single',...
%     'Name','Select source cluster',...
%     'ListString',groupnames);
% if isempty(s1) || ~v1
%     disp('*** WARNING: Viewing loglikelihood cancelled.');
%     return
% elseif s1==npops+1
%     view_all_density(changesInLogml, partition, filename);
% else
%     remain_pop = logical((1:npops)~=s1);
%     %waitALittle;
%     [s2,v2] = listdlg('PromptString', 'Select target clusters:',...
%         'SelectionMode','multiple',...
%         'Name','Select target cluster',...
%         'ListString',groupnames(remain_pop));
%     if isempty(s2) || ~v2
%         disp('*** WARNING: Viewing loglikelihood cancelled.');
%         return
%     end
%     inds = logical(partition==s1);
%     ninds = sum(inds); % individuals in the source cluster
%     remain = find(remain_pop);
%     fprintf('Source cluster: %d\n', s1);
%     fprintf('Number of strains: %d\n',ninds);
%     fprintf('Target cluster(s): %s\n', num2str(remain(s2)));
% 
%     f = zeros(length(s2),100);
%     xi = zeros(length(s2),100);
%     for i = 1:length(s2)
%         [f(i,:),xi(i,:)] = ksdensity_myown(changesInLogml(remain(s2(i)),inds)');
%     end
% 
%     map = giveColors(npops);
%     h0 = figure('NumberTitle', 'off'); % density plot;
%     set(h0,'Tag','density_plot');
%     set(h0,'Name',['Density of log likelihood changes - ' filename]);
%     set(gca,'ColorOrder',map(remain(s2),:));
%     hold on;
%     plot(xi',f');
%     set(gca,'xlim',[min(min(xi)),max(max(xi))]);
%     legend(groupnames(remain(s2)));
%     xlabel('Change of log likelihood');
%     ylabel('Estimated density');
%     title(sprintf('Cluster %d',s1));
% 
%     %     h1 = figure('NumberTitle', 'off'); % histogram plot;
%     %     set(h1,'Tag','histogram_plot');
%     %     set(h1,'Name',['Histogram of log likelihood changes - ' filename]);
%     %     hist(changesInLogml(remain(s2),inds)');
%     %     colormap(map(remain(s2),:));
%     %     % histfit(changesInLogml(remain(s2),inds)'); % need statistical toolbox
%     %     legend(groupnames(remain(s2)));
%     %     xlabel('Change of log likelihood');
%     %     ylabel('Frequency');
%     %     title(sprintf('Cluster %d',s1));
%     % rose(changesInLogml(remain(s2),inds));
% end


%--------------------------------------%
%%% functions for linkage model %%%
%---------------------------------------

%--------------------------------------------------------------------------

function changes = computeChanges(cq_counts, cq_sumcounts,...
    sp_counts, sp_sumcounts,...
    partition, pop_logml,...
    ind, adjprior_cq, adjprior_sp, ...
    indCqCounts, indSpCounts)
% Computes changes in log-marginal likelihood if individual ind is
% moved to another population
%
% Input:
% ind - the individual to be moved
% adjprior_cq & _sp  - adjpriors for cliques and separators
% indCqCounts, indSpCounts - counts for individual ind
%
% Output:
% changes - table of size 1*npops. changes(i) = difference in logml if
% ind is move to population i.

npops = size(cq_counts,3);
changes = zeros(npops,1);

i1 = partition(ind);
i1_logml = pop_logml(i1);
sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

cq_counts(:,:,i1) = cq_counts(:,:,i1)-indCqCounts;
cq_sumcounts(i1,:) = cq_sumcounts(i1,:)-sumCq;
sp_counts(:,:,i1) = sp_counts(:,:,i1)-indSpCounts;
sp_sumcounts(i1,:) = sp_sumcounts(i1,:)-sumSp;

new_i1_logml = computePopulationLogml(double(cq_counts), double(cq_sumcounts),...
    double(sp_counts), double(sp_sumcounts),...
    i1, adjprior_cq, adjprior_sp);

cq_counts(:,:,i1) = cq_counts(:,:,i1)+indCqCounts;
cq_sumcounts(i1,:) = cq_sumcounts(i1,:)+sumCq;
sp_counts(:,:,i1) = sp_counts(:,:,i1)+indSpCounts;
sp_sumcounts(i1,:) = sp_sumcounts(i1,:)+sumSp;


i2 = [1:i1-1 , i1+1:npops];
i2_logml = pop_logml(i2);

cq_counts(:,:,i2) = cq_counts(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
cq_sumcounts(i2,:) = cq_sumcounts(i2,:)+repmat(sumCq,[npops-1 1]);
sp_counts(:,:,i2) = sp_counts(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
sp_sumcounts(i2,:) = sp_sumcounts(i2,:) + repmat(sumSp,[npops-1 1]);

new_i2_logml = computePopulationLogml(double(cq_counts), double(cq_sumcounts),...
    double(sp_counts), double(sp_sumcounts),...
    i2, adjprior_cq, adjprior_sp);

changes(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%--------------------------------------------------------------------------

function popLogml = computePopulationLogml(cq_counts, cq_sumcounts,...
    sp_counts, sp_sumcounts,...
    pops, adjprior_cq, adjprior_sp)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on määritelty pops-muuttujalla.


nall_cq = size(cq_counts,1);
nall_sp = size(sp_counts, 1);
ncliq = size(cq_counts,2);
nsep = size(sp_counts, 2);

z = length(pops);

popLogml_cq = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_cq,[1 1 length(pops)]) + cq_counts(:,:,pops)) ...
    ,[nall_cq ncliq z]),1),2)) - sum(gammaln(1+cq_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_cq)));

popLogml_sp = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_sp,[1 1 length(pops)]) + sp_counts(:,:,pops)) ...
    ,[nall_sp nsep z]),1),2)) - sum(gammaln(1+sp_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_sp)));

popLogml = popLogml_cq - popLogml_sp;
clear cq_counts cq_sumcounts sp_counts sp_sumcounts;

%--------------------------------------------------------------------------


%--------------------------------------%
%%% functions for independence model %%%
%---------------------------------------
%--------------------------------------------------------------------------
function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, rowsFromInd, data, adjprior, priorTerm, ...
                   COUNTS, SUMCOUNTS, PARTITION, POP_LOGML)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
% muutos logml:ss? mikäli yksil?ind siirretään koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lisättäv?
% COUNTS:in siivuun i2, mikäli muutos toteutetaan.

npops = size(COUNTS,3);
muutokset = zeros(npops,1);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);

rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml2(i1, adjprior, priorTerm,COUNTS, SUMCOUNTS);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml2(i2, adjprior, priorTerm,COUNTS, SUMCOUNTS);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;

%--------------------------------------------------------------------------
function popLogml = computePopulationLogml2(pops, adjprior, priorTerm, ...
                                            COUNTS, SUMCOUNTS)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on määritelty pops-muuttujalla.

x = size(COUNTS,1);
y = size(COUNTS,2);
z = length(pops);

popLogml = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
    ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;


%--------------------------------------------------------------------------
function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukumäärät (vastaavasti kuin COUNTS:issa), jotka ovat data:n 
% riveill?rows.

diffInCounts = zeros(max_noalle, nloci);
for i=rows
    row = data(i,:);
    notEmpty = find(row>=0);
    
    if length(notEmpty)>0
        diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
            diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
    end
end    


