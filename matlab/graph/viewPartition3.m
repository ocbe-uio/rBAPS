function [partition,data,popnames_return,pvalue] = ...
           viewPartition3(data, popnames, npops, admixnpops, ...
                          inliers, partition, pvalue, filename)
% Kuin viewParition.m, mutta toimii itsenäisesti. Ei siis tarvitse
% globaalia COUNTS-muuttujaa.

%global COUNTS;
%npops = size(COUNTS, 3);

% Determine the order of individuals by calculating the distances.
nind = size(data,1);
% inliers = unique(partition);
% inliers = find(any(data));

dist = zeros(1,nind*(nind-1)/2);
k = 1;
% m = 1;
disp('---------------------------------------------------');
fprintf('Permutating the individuals according to their relateness.\n');
fprintf('This may take several minutes. Please wait...');
for i=1:nind-1
    for j=i+1:nind
        dist(k) = distance(data(i,:)',data(j,:)');
        k = k +1;
    end
end
Z = linkage(dist);
[T,perm] = dendrogram_alpha(Z,0);
% [H,T,perm] = dendrogram_beta(Z,popnames(:,1)); % This function gives dendrograms.
data = data(perm,:);
popnames = popnames(perm,:);
popnames_return = popnames; % this is the popnames to be saved
partition = partition(perm,:);
pvalue = pvalue(perm,:);

% npairs = nchoosek(nind,2);
% pairwisesimilarity = zeros(npairs,1);
% clusters = cell(1,ncluster);
% for i = 1:ncluster
%     clusters{i} = find(data(:,i)>0);
% end
% % 'adjacency matrix'
% adj_mat = cell(ncluster,ncluster);
% k = 2;
% for i = 1:ncluster
%     for j = k:ncluster
%         adj_mat{i,j} = intersect(clusters{i},clusters{j});
%     end
%     k = k + 1;
% end
% 
% adj_mat = zeros(ncluster,ncluster);
% k = 2;
% for i = 1:ncluster
%     for j = k:ncluster
%         adj_mat(i,j) = length(intersect(clusters{i},clusters{j}));
%     end
%     k = k + 1;
% end
% adj_mat = adj_mat+adj_mat';
% admix_count = sum(adj_mat,1);
% order = zeros(1,ncluster);
% middle_cluster = find(recomb_count==max(admix_count));
% reference_cluster = middle_cluster;
% middle_point = floor(ncluster/2);
% order(middle_point) = middle_cluster(1);
% for i = middle_point+1:ncluster
%     order(i) = find(adj_mat(:,reference_cluster)==max(adj_mat(:,reference_cluster)));
%     adj_mat(referecen_cluster,order(i)) = 0;
%     adj_mat(order(i),reference_cluster) = 0;
%     reference_cluster = order(i);
% end
% [ordered_cluster,ix] = sort(size_cluster,'descend');

% disp(['Number of populations: ' num2str(npops)]);
% if npops>30
%     disp(' ');
%     disp('Figure can be drawn only if the number of populations');
%     disp('is less or equal to 30.');
%     disp(' ');
%     return;
% end


varit = giveColors(npops);
korkeinviiva = 1.05;
pieninarvo = -korkeinviiva;

h0 = figure;
set(h0, 'NumberTitle', 'off'); %image_figure;   %Muutettu
set(h0,'Tag','admixture_result');
set(h0,'Name',['Admixture Result After Permutation - ' filename]); 

if size(popnames,1) > 150 % display only groups
    for i=1:nind
        popnames{i,1} = cell({'|'});  % CELL in the CELL
    end
end

for i = inliers
    pool = find(partition==i);
    if isempty(pool)
        continue
    end
    chosen_ix = 1 + floor(length(pool)/2);
%     if ~find(pool==chosen_one) % gap, find the other way
%         chosen_one = pool(end) - floor(length(pool)/2);
%         if ~find(pool==chosen_one) % then the algorithm really sucks
%            chosen_one = pool(1);
%         end
%     end
    chosen_one = pool(chosen_ix);
    popnames{chosen_one,1} = {[popnames{chosen_one,1}{1} sprintf('Cluster%d',i)]};
end
   
   
tiedot.popnames = popnames;
tiedot.info = data;
set(h0,'UserData',tiedot);

set(gca,'Xlim', [-.5 ,nind+.5], 'YLim', [pieninarvo ,korkeinviiva], ...
    'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

% Display the color bars
for i=1:nind
    
    if any(data(i,:)>0)
        cumOsuudet = cumsum(data(i,:));
    
        % Pylvään piirtäminen
        for j=1:admixnpops
            if j==1
                if cumOsuudet(1)>0
                    h0 =patch([i-1, i, i, i-1], [0, 0, cumOsuudet(1), cumOsuudet(1)], varit(inliers(j),:));
                    set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
                end
            else
                if (cumOsuudet(j)>cumOsuudet(j-1))
                    h0 = patch([i-1, i, i, i-1], [cumOsuudet(j-1), cumOsuudet(j-1), ...
                    cumOsuudet(j), cumOsuudet(j)], varit(inliers(j),:));
                    set(h0,'EdgeColor','none'); % Midevaa varten kommentoitava!
                end
            end
        end
    end
end


% Display the texts
if ~isempty(popnames)
    npops = size(popnames,1); % NB! npops now stands for the number of individuals.
    for i=1:npops
        firstInd = popnames{i,2};
        if size(popnames,1)~=nind
            line([firstInd-1, firstInd-1], [0,1], 'Color', 'k');  %Populaatioiden rajat
        end
        
        % The determination of x_paikka is changed, since popnames could be
        % non-sequential. - Jing
        if i<npops
            % x_paikka = popnames{i,2}-1+(popnames{i+1,2}-popnames{i,2})/2;
            x_paikka = i-1;
        else
            % x_paikka = popnames{i,2}-1+(nind+1-popnames{i,2})/2;
            x_paikka = i-1+(nind+1-i)/2;
        end
               
        korkeuskerroin = pieninarvo / -0.2;
        suhdekerroin = npops/6;
        names = popnames{i,1}{1};
        max_length = length(names);
        for letter_num = 1: max_length % {1}
            letter = names(letter_num);%alter .004|
            color_ix = find(data(i,:)==max(data(i,:)));
            color_ix = color_ix(1);
            text(x_paikka+korjaus(letter)*suhdekerroin, ...
                0.0005*korkeuskerroin-0.02*letter_num*korkeuskerroin, ...
                letter, 'Interpreter','none','Color',varit(inliers(color_ix),:));
%          text(x_paikka+korjaus(letter)*suhdekerroin, ...
%                 pieninarvo*letter_num/max_length, ...
%                 letter, 'Interpreter','none','Color',varit(color_ix,:));
        end
    end
    line([nind,nind],[0,1],'Color','k');
end
fprintf('Finished.\n');
%-------------------------------------------------------------------------------------

function extra = korjaus(letter)
if any(letter == 'ijlI')
    extra = 0.022;
elseif any(letter == 'r')
    extra = 0.016;
elseif any(letter == 'k')
    extra = 0.009;
elseif any(letter == 'f')
    extra = 0.013;
elseif any(letter == 't')
    extra = 0.014;
elseif any(letter == 'w')
    extra = -0.003;
else
    extra = 0;
end;