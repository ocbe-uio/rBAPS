function mutationPlot(tietue)

global COUNTS; global PARTITION; global SUMCOUNTS;
clearGlobalVars;

h0 = findobj('Tag','load_menu');
c = get(h0,'UserData');
h0 = findobj('Tag','filename1_text');
filename1 = get(h0,'String');
disp('---------------------------------------------------');
disp('Viewing the mutation plot.');
disp(['Load the mixture result from: ',[filename1],'...']);

if (~isstruct(tietue))
%     [filename, pathname] = uigetfile('*.mat', 'Load mixture result file');
%     if (filename==0 & pathname==0) return; 
%     else
%         %display('---------------------------------------------------');
%         %display(['Reading mixture result from: ',[pathname filename],'...']);
%     end
%     pause(0.0001);
%     h0 = findobj('Tag','filename1_text');
%     set(h0,'String',filename); clear h0;
%     
%     struct_array = load([pathname filename]);
%     if isfield(struct_array,'c')  %Matlab versio
%         c = struct_array.c;
%         if ~isfield(c,'PARTITION') | ~isfield(c,'rowsFromInd')
%             disp('Incorrect file format');
%             return
%         end
%     else
%         disp('Incorrect file format');
%         return;
%     end
    PARTITION = c.PARTITION; COUNTS = c.COUNTS; SUMCOUNTS = c.SUMCOUNTS;
    alleleCodes = c.alleleCodes; adjprior = c.adjprior; popnames = c.popnames;
    rowsFromInd = c.rowsFromInd; data = c.data; npops = c.npops; noalle = c.noalle;
else
    PARTITION = tietue.PARTITION;
    COUNTS = tietue.COUNTS;
    SUMCOUNTS = tietue.SUMCOUNTS;
    alleleCodes = tietue.alleleCodes;
    adjprior = tietue.adjprior;
    popnames = tietue.popnames;
    rowsFromInd = tietue.rowsFromInd;
    data = double(tietue.data);
    npops = tietue.npops;
    noalle = tietue.noalle;
end

nloci = size(COUNTS,2);
ninds = size(data,1)/rowsFromInd;

if isfield(c,'CQ_COUNTS')
    % Linked data
    if isfield(c,'gene_lengths')
        gene_lengths = c.gene_lengths;
    else
        [filename, pathname] = uigetfile('*.txt', 'Load file with gene lengths.');
        if (filename==0 & pathname==0) 
            return 
        end
        gene_lengths = load([pathname filename]);
    end
    
    component_mat = zeros(length(gene_lengths), max(gene_lengths));
    cum_length = cumsum(gene_lengths);
    component_mat(1,1:gene_lengths(1))=1:gene_lengths(1);
    for i = 2:length(gene_lengths)
        component_mat(i,1:gene_lengths(i)) = cum_length(i-1)+1:cum_length(i);
    end
else
    component_mat = 1:nloci;
    gene_lengths = nloci;
end

if isfield(c,'CQ_COUNTS')
    answers = inputdlg({'Index of the individual of interest',...
        'BF limit (log)',...
        'Genes to analyze'},...
        'INPUT',[1; 1; 1],...
        {' ','2.3',['1:' num2str(length(gene_lengths))]});
    ind = str2num(answers{1});
    BF = str2num(answers{2});
    genes = str2num(answers{3});
    %ind = input('Input individual: ');
    %BF = input('Input BF: ');
    %genes = input('Input genes to analyze: ');
    all_right = check_inputs(ind,BF,genes, ninds, gene_lengths);
    if all_right==0
        return
    end
    n_genes = length(genes);
    nameText = ['ind: ' num2str(ind) ', genes: ' num2str(genes)];
else
    answers = inputdlg({'Index of the individual of interest',...
        'BF limit (log)'},'INPUT',[1; 1],{' ','2.3'});
    ind = str2num(answers{1});
    BF = str2num(answers{2});
    %ind = input('Input individual: ');
    %BF = input('Input BF: ');
    all_right = check_inputs(ind,BF,1, ninds, gene_lengths);
    if all_right==0
        return
    end
    genes = 1;
    n_genes = 1;
    nameText = ['ind: ' num2str(ind)];
end

origin = PARTITION(ind);
varit = giveColors(npops);

figure('NumberTitle','off','Name',nameText);

mutationList = repmat(...
    struct('gene',[],'site',[],...
    'row',[],'pops',[],'bf_vals',[]), [10 1]);
n_mutations = 0;

for i = 1:n_genes
    % Loop over different genes
    gene = genes(i);
    
    for j = 1:rowsFromInd
        % Loop over different "haplotypes"
        left = 0.05;
        bottom = 1 - (1/n_genes)*i + (1/n_genes)/(rowsFromInd*2+1)*(2*j-1);
        width = 0.9;
        height = 1/n_genes/(rowsFromInd*2+1);
        axes('Position',[left bottom width height]);
        set(gca, 'Xlim', [-.5 , gene_lengths(gene)+.5], 'YLim', [0,1], ...
            'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
        
        for k=1:gene_lengths(gene)
            % Loop over different sites
            site = component_mat(gene,k);
            allele = data((ind-1)*rowsFromInd+j, site);
            %disp(num2str((ind-1)*rowsFromInd+j));
            counts = squeeze(COUNTS(:,site,:));
            
            logml = zeros(npops,1);
            if allele>0
                for m=1:npops
                    % Calculate the predictive probabilities of different
                    % origins
                    counts(allele,origin) = counts(allele,origin)-1;
                    counts(allele,m) = counts(allele,m)+1;
                    logml(m) = calculateLogml(counts, noalle(site));
                    counts(allele,m) = counts(allele,m)-1;
                    counts(allele,origin) = counts(allele,origin)+1;
                end
            end
            
            logml = logml-max(logml);
            probs = exp(logml) ./ sum(exp(logml));
            
            cumprobs = cumsum(probs);
            
            if log(max(probs)/probs(origin)) >= BF
                h0=patch([k-1 k k k-1], [0 0 cumprobs(1) cumprobs(1)], varit(1,:));
                set(h0,'EdgeColor','none');
                for m=2:npops
                    h0=patch([k-1 k k k-1], [cumprobs(m-1) cumprobs(m-1) cumprobs(m) cumprobs(m)], varit(m,:));
                    set(h0,'EdgeColor','none');
                end
                n_mutations = n_mutations+1;
                if n_mutations>length(mutationList)
                    mutationList = [mutationList; repmat(...
                        struct('gene',[],'site',[],...
                        'row',[],'pops',[],'bf_vals',[]), [length(mutationList) 1])];
                end
                mutationList(n_mutations).gene = gene;
                mutationList(n_mutations).site = k; % The location in gene i!
                mutationList(n_mutations).row = j;
                aux = log(probs ./ probs(origin));
                mutationList(n_mutations).pops = find(aux>=BF);
                mutationList(n_mutations).bf_vals = aux(find(aux>=BF));
            end
        end
    end
end
mlst_data = isfield(c,'CQ_COUNTS');
writeResults(mutationList, n_mutations, ind, origin, mlst_data);


%--------------------------------------------------------------------


function all_right = check_inputs(ind,BF,genes, ninds, gene_lengths)

all_right = 0;
if length(ind)~=1
    disp('ERROR: index of one individual must be given.');
    return;
end
if ind<=0 | ind >ninds
    disp('ERROR: Index of the given individual is out of range.');
    return
end
if length(BF)~=1
    disp('ERROR: one BF value must be given.');
    return
end
if BF<0
    disp('ERROR: BF must be positive.');
    return
end
if length(genes)<1 | length(genes)>length(gene_lengths)
    disp('ERROR: input for the genes was incorrect.');
    return
end
if any(genes<1) | any(genes>length(gene_lengths))
    disp('ERROR: input for the genes was incorrect.');
    return
end
all_right = 1;

%------------------------------------------------------------------

function writeResults(mutationList, n_mutations, ind, home, mlst_data)

h0 = findobj('Tag','filename2_text');
outf = get(h0,'String'); clear h0;

if length(outf)>0
    fid = fopen(outf,'a');
else
    fid = -1;
end

disp(' ');
disp('--------------------------');
disp(['Individual: ' num2str(ind) '.']);
disp(['Origin of the individual ' num2str(home) '.']);

if fid ~= -1
    fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['--------------------------------------------']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Individual: ' num2str(ind) '.']); fprintf(fid, '\n');
    fprintf(fid,'%s \n', ['Origin of the individual ' num2str(home) '.']); fprintf(fid, '\n');
end


if n_mutations==0
    disp('No putative locations for mutations detected');
    if fid ~= -1
        fprintf(fid,'%s \n', 'No putative locations for mutations detected'); fprintf(fid, '\n');
    end
    return
end

disp(' ');
if fid ~= -1
    fprintf(fid,'%s \n', ['  ']); fprintf(fid, '\n');
end

if mlst_data==1
    disp('Putative locations for mutations:');
    disp('gene,  site,   possible origins');
    if fid ~= -1
        fprintf(fid,'%s \n', 'Putative locations for mutations:'); fprintf(fid, '\n');
        fprintf(fid,'%s \n', 'gene,  site,   possible origins'); fprintf(fid, '\n');
    end
    for i=1:n_mutations
        gene = mutationList(i).gene;
        site = mutationList(i).site;
        text_line = [num2str(gene) blanks(6-floor(log10(gene))) num2str(site) blanks(7-floor(log10(site)))];
        pops = mutationList(i).pops;
        bf_vals = mutationList(i).bf_vals;
        for j=1:length(pops)
            text_line = [text_line num2str(pops(j)) '(' num2str(bf_vals(j)) ') '];
        end
        disp(text_line);
        if fid ~= -1
            fprintf(fid,'%s \n', [text_line]); fprintf(fid, '\n');
        end
    end
else
    disp('Putative locations for mutations:');
    disp('locus,  haplotype, possible origins');
    if fid ~= -1
        fprintf(fid,'%s \n', 'Putative locations for mutations:'); fprintf(fid, '\n');
        fprintf(fid,'%s \n', 'locus,  haplotype, possible origins'); fprintf(fid, '\n');
    end
    for i=1:n_mutations
        locus = mutationList(i).site;
        row = mutationList(i).row;
        text_line = [num2str(locus) blanks(7-floor(log10(locus))) num2str(row) blanks(10-floor(log10(row)))];
        pops = mutationList(i).pops;
        bf_vals = mutationList(i).bf_vals;
        for j=1:length(pops)
            text_line = [text_line num2str(pops(j)) '(' num2str(bf_vals(j)) ') '];
        end
        disp(text_line);
        if fid ~= -1
            fprintf(fid,'%s \n', [text_line]); fprintf(fid, '\n');
        end
    end
end

if fid ~= -1
    fclose(fid);
end

%------------------------------------------------------------------

function val = calculateLogml(counts, noalle)
% counts corresponds to counts of a SINGLE locus. It is two-dimensional
% with as many columns as there are populations.

npops = size(counts,2);

prior = ones(noalle,npops)./noalle;
counts = counts(1:noalle,:);

val = sum(gammaln(sum(prior,1)),2) ...
    - sum(gammaln(sum(counts+prior,1)),2) ...
    + sum(sum(gammaln(counts+prior))) ...
    - sum(sum(gammaln(prior)));


%---------------------------------------------------------------

function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];