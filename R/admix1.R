#' @title Admixture analysis
#' @param tietue record
#' @details If the record == -1, the mixture results file is loaded. Otherwise, will the required variables be retrieved from the record fields?
#' @export
admix1 <- function(tietue) {
    if (!is.list(tietue)) {
#         c(filename, pathname) = uigetfile('*.mat', 'Load mixture result file');
#         if (filename==0 & pathname==0), return; 
#         else
#             disp('---------------------------------------------------');
#             disp(['Reading mixture result from: ',[pathname filename],'...']);
#         end
#         pause(0.0001);
#         h0 = findobj('Tag','filename1_text');
#         set(h0,'String',filename); clear h0;
        
#         struct_array = load([pathname filename]);
#         if isfield(struct_array,'c')  #Matlab versio
#             c = struct_array.c;
#             if ~isfield(c,'PARTITION') | ~isfield(c,'rowsFromInd')
#                 disp('Incorrect file format');
#                 return
#             end
#         elseif isfield(struct_array,'PARTITION')  #Mideva versio
#             c = struct_array;
#             if ~isfield(c,'rowsFromInd')
#                 disp('Incorrect file format');
#                 return
#             end
#         else
#             disp('Incorrect file format');
#             return;
#         end
        
#         if isfield(c, 'gene_lengths') && ...
#                 (strcmp(c.mixtureType,'linear_mix') | ...
#                 strcmp(c.mixtureType,'codon_mix')) # if the mixture is from a linkage model
#             # Redirect the call to the linkage admixture function.
#             c.data = noIndex(c.data,c.noalle); # call function noindex to remove the index column
#             linkage_admix(c);
#             return
#         end
        
#         PARTITION = c.PARTITION; COUNTS = c.COUNTS; SUMCOUNTS = c.SUMCOUNTS;
#         alleleCodes = c.alleleCodes; adjprior = c.adjprior; popnames = c.popnames;
#         rowsFromInd = c.rowsFromInd; data = c.data; npops = c.npops; noalle = c.noalle;
    } else {
        PARTITION <- tietue$PARTITION
        COUNTS <- tietue$COUNTS
        SUMCOUNTS <- tietue$SUMCOUNTS
        alleleCodes <- tietue$alleleCodes
        adjprior <- tietue$adjprior
        popnames <- tietue$popnames
        rowsFromInd <- tietue$rowsFromInd
        data <- as.double(tietue$data)
        npops <- tietue$npops
        noalle <- tietue$noalle
    }

#     answers = inputdlg({['Input the minimum size of a population that will'...
#                 ' be taken into account when admixture is estimated.']},...
#                 'Input minimum population size',[1],...
#                 {'5'});
#     if isempty(answers) return; end
#     alaRaja = str2num(answers{1,1});
#     [npops] = poistaLiianPienet(npops, rowsFromInd, alaRaja);

#     nloci = size(COUNTS,2);
#     ninds = size(data,1)/rowsFromInd;

#     answers = inputdlg({['Input number of iterations']},'Input',[1],{'50'});
#     if isempty(answers) return; end
#     iterationCount = str2num(answers{1,1});

#     answers = inputdlg({['Input number of reference individuals from each population']},'Input',[1],{'50'});
#     if isempty(answers) nrefIndsInPop = 50;
#     else nrefIndsInPop = str2num(answers{1,1});
#     end

#     answers = inputdlg({['Input number of iterations for reference individuals']},'Input',[1],{'10'});
#     if isempty(answers) return; end
#     iterationCountRef = str2num(answers{1,1});


#     # First calculate log-likelihood ratio for all individuals:
#     likelihood = zeros(ninds,1);
#     allfreqs = computeAllFreqs2(noalle);
#     for ind = 1:ninds
#         omaFreqs = computePersonalAllFreqs(ind, data, allfreqs, rowsFromInd);
#         osuusTaulu = zeros(1,npops);
#         if PARTITION(ind)==0
#             # Yksil?on outlier
#         elseif PARTITION(ind)~=0
#             if PARTITION(ind)>0
#                 osuusTaulu(PARTITION(ind)) = 1;
#             else
#                 # Yksilöt, joita ei ole sijoitettu mihinkään koriin.
#                 arvot = zeros(1,npops);
#                 for q=1:npops
#                     osuusTaulu = zeros(1,npops);
#                     osuusTaulu(q) = 1;
#                     arvot(q) = computeIndLogml(omaFreqs, osuusTaulu);
#                 end
#                 [iso_arvo, isoimman_indeksi] = max(arvot);
#                 osuusTaulu = zeros(1,npops);
#                 osuusTaulu(isoimman_indeksi) = 1;
#                 PARTITION(ind)=isoimman_indeksi;
#             end
#             logml = computeIndLogml(omaFreqs, osuusTaulu);
#             logmlAlku = logml;
#             for osuus = [0.5 0.25 0.05 0.01]
#                 [osuusTaulu, logml] = etsiParas(osuus, osuusTaulu, omaFreqs, logml);
#             end
#             logmlLoppu = logml;
#             likelihood(ind) = logmlLoppu-logmlAlku;
#         end
#     end

#     # Analyze further only individuals who have log-likelihood ratio larger than 3:
#     to_investigate = (find(likelihood>3))';
#     disp('Possibly admixed individuals: ');
#     for i = 1:length(to_investigate)
#         disp(num2str(to_investigate(i)));
#     end
#     disp(' ');
#     disp('Populations for possibly admixed individuals: ');
#     admix_populaatiot = unique(PARTITION(to_investigate));
#     for i = 1:length(admix_populaatiot)
#         disp(num2str(admix_populaatiot(i)));
#     end

#     # THUS, there are two types of individuals, who will not be analyzed with
#     # simulated allele frequencies: those who belonged to a mini-population
#     # which was removed, and those who have log-likelihood ratio less than 3.
#     # The value in the PARTITION for the first kind of individuals is 0. The
#     # second kind of individuals can be identified, because they do not
#     # belong to "to_investigate" array. When the results are presented, the
#     # first kind of individuals are omitted completely, while the second kind
#     # of individuals are completely put to the population, where they ended up
#     # in the mixture analysis. These second type of individuals will have a
#     # unit p-value.


#     # Simulate allele frequencies a given number of times and save the average
#     # result to "proportionsIt" array.

#     proportionsIt = zeros(ninds,npops);
#     for iterationNum = 1:iterationCount
#         disp(['Iter: ' num2str(iterationNum)]);
#         allfreqs = simulateAllFreqs(noalle);   # Allele frequencies on this iteration.
        
#         for ind=to_investigate
#             #disp(num2str(ind));
#             omaFreqs = computePersonalAllFreqs(ind, data, allfreqs, rowsFromInd);
#             osuusTaulu = zeros(1,npops);
#             if PARTITION(ind)==0
#                 # Yksil?on outlier
#             elseif PARTITION(ind)~=0
#                 if PARTITION(ind)>0
#                     osuusTaulu(PARTITION(ind)) = 1;
#                 else
#                     # Yksilöt, joita ei ole sijoitettu mihinkään koriin.
#                     arvot = zeros(1,npops);
#                     for q=1:npops
#                         osuusTaulu = zeros(1,npops);
#                         osuusTaulu(q) = 1;
#                         arvot(q) = computeIndLogml(omaFreqs, osuusTaulu);
#                     end
#                     [iso_arvo, isoimman_indeksi] = max(arvot);
#                     osuusTaulu = zeros(1,npops);
#                     osuusTaulu(isoimman_indeksi) = 1;
#                     PARTITION(ind)=isoimman_indeksi;
#                 end
#                 logml = computeIndLogml(omaFreqs, osuusTaulu);
                
#                 for osuus = [0.5 0.25 0.05 0.01]
#                     [osuusTaulu, logml] = etsiParas(osuus, osuusTaulu, omaFreqs, logml);
#                 end
#             end
#             proportionsIt(ind,:) = proportionsIt(ind,:).*(iterationNum-1) + osuusTaulu;
#             proportionsIt(ind,:) = proportionsIt(ind,:)./iterationNum;
#         end
#     end

#     #disp(['Creating ' num2str(nrefIndsInPop) ' reference individuals from ']);
#     #disp('each population.');

#     #allfreqs = simulateAllFreqs(noalle);  # Simuloidaan alleelifrekvenssisetti
#     allfreqs = computeAllFreqs2(noalle); # Koitetaan tällaista.


#     # Initialize the data structures, which are required in taking the missing
#     # data into account:
#     n_missing_levels = zeros(npops,1);  # number of different levels of "missingness" in each pop (max 3).
#     missing_levels = zeros(npops,3);    # the mean values for different levels.
#     missing_level_partition = zeros(ninds,1); # level of each individual (one of the levels of its population).
#     for i=1:npops
#         inds = find(PARTITION==i);
#         # Proportions of non-missing data for the individuals:
#         non_missing_data = zeros(length(inds),1);
#         for j = 1:length(inds)
#             ind = inds(j);
#             non_missing_data(j) = length(find(data((ind-1)*rowsFromInd+1:ind*rowsFromInd,:)>0)) ./ (rowsFromInd*nloci);
#         end
#         if all(non_missing_data>0.9)
#             n_missing_levels(i) = 1;
#             missing_levels(i,1) = mean(non_missing_data);
#             missing_level_partition(inds) = 1;
#         else
#             [ordered, ordering] = sort(non_missing_data);
#             #part = learn_simple_partition(ordered, 0.05);
#             part = learn_partition_modified(ordered);
#             aux = sortrows([part ordering],2);
#             part = aux(:,1);
#             missing_level_partition(inds) = part;
#             n_levels = length(unique(part));
#             n_missing_levels(i) = n_levels;
#             for j=1:n_levels
#                 missing_levels(i,j) = mean(non_missing_data(find(part==j)));
#             end
#         end
#     end

#     # Create and analyse reference individuals for populations
#     # with potentially admixed individuals:
#     refTaulu = zeros(npops,100,3);
#     for pop = admix_populaatiot'

#         for level = 1:n_missing_levels(pop)
            
#             potential_inds_in_this_pop_and_level = ...
#                 find(PARTITION==pop & missing_level_partition==level &...
#                 likelihood>3);  # Potential admix individuals here.
            
#             if ~isempty(potential_inds_in_this_pop_and_level)
            
#                 #refData = simulateIndividuals(nrefIndsInPop,rowsFromInd,allfreqs);
#                 refData = simulateIndividuals(nrefIndsInPop, rowsFromInd, allfreqs, ...
#                     pop, missing_levels(pop,level));
            
#                 disp(['Analysing the reference individuals from pop ' num2str(pop) ' (level ' num2str(level) ').']);
#                 refProportions = zeros(nrefIndsInPop,npops);
#                 for iter = 1:iterationCountRef
#                     #disp(['Iter: ' num2str(iter)]);
#                     allfreqs = simulateAllFreqs(noalle);
                
#                     for ind = 1:nrefIndsInPop
#                         omaFreqs = computePersonalAllFreqs(ind, refData, allfreqs, rowsFromInd);
#                         osuusTaulu = zeros(1,npops);
#                         osuusTaulu(pop)=1;
#                         logml = computeIndLogml(omaFreqs, osuusTaulu);
#                         for osuus = [0.5 0.25 0.05 0.01]
#                             [osuusTaulu, logml] = etsiParas(osuus, osuusTaulu, omaFreqs, logml);
#                         end
#                         refProportions(ind,:) = refProportions(ind,:).*(iter-1) + osuusTaulu;
#                         refProportions(ind,:) = refProportions(ind,:)./iter;
#                     end
#                 end
#                 for ind = 1:nrefIndsInPop
#                     omanOsuus = refProportions(ind,pop);
#                     if round(omanOsuus*100)==0
#                         omanOsuus = 0.01;
#                     end
#                     if abs(omanOsuus)<1e-5
#                         omanOsuus = 0.01;
#                     end
#                     refTaulu(pop, round(omanOsuus*100),level) = refTaulu(pop, round(omanOsuus*100),level)+1;
#                 end
#             end
#         end
#     end

#     # Rounding of the results:
#     proportionsIt = proportionsIt.*100; proportionsIt = round(proportionsIt);
#     proportionsIt = proportionsIt./100;
#     for ind = 1:ninds
#         if ~any(to_investigate==ind)
#             if PARTITION(ind)>0
#                 proportionsIt(ind,PARTITION(ind))=1;
#             end
#         else
#             # In case of a rounding error, the sum is made equal to unity by
#             # fixing the largest value.
#             if (PARTITION(ind)>0) & (sum(proportionsIt(ind,:)) ~= 1)
#                 [isoin,indeksi] = max(proportionsIt(ind,:));
#                 erotus = sum(proportionsIt(ind,:))-1;
#                 proportionsIt(ind,indeksi) = isoin-erotus;
#             end
#         end
#     end

#     # Calculate p-value for each individual:
#     uskottavuus = zeros(ninds,1);
#     for ind = 1:ninds
#         pop = PARTITION(ind);
#         if pop==0 # Individual is outlier
#             uskottavuus(ind)=1;
#         elseif isempty(find(to_investigate==ind))
#             # Individual had log-likelihood ratio<3
#             uskottavuus(ind)=1;
#         else
#             omanOsuus = proportionsIt(ind,pop);
#             if abs(omanOsuus)<1e-5
#                 omanOsuus = 0.01;
#             end
#             if round(omanOsuus*100)==0
#                 omanOsuus = 0.01;
#             end
#             level = missing_level_partition(ind);
#             refPienempia = sum(refTaulu(pop, 1:round(100*omanOsuus), level));
#             uskottavuus(ind) = refPienempia / nrefIndsInPop;
#         end
#     end

#     tulostaAdmixtureTiedot(proportionsIt, uskottavuus, alaRaja, iterationCount); 

#     viewPartition(proportionsIt, popnames);

#     talle = questdlg(['Do you want to save the admixture results?'], ...
#         'Save results?','Yes','No','Yes');
#     if isequal(talle,'Yes')
#         #waitALittle;
#         [filename, pathname] = uiputfile('*.mat','Save results as');


#         if (filename == 0) & (pathname == 0)
#             # Cancel was pressed
#             return
#         else # copy 'baps4_output.baps' into the text file with the same name.
#             if exist('baps4_output.baps','file')
#                 copyfile('baps4_output.baps',[pathname filename '.txt'])
#                 delete('baps4_output.baps')
#             end
#         end


#         if (~isstruct(tietue))
#             c.proportionsIt = proportionsIt; 
#             c.pvalue = uskottavuus; # Added by Jing
#             c.mixtureType = 'admix'; # Jing
#             c.admixnpops = npops;
#     #         save([pathname filename], 'c');
#             save([pathname filename], 'c', '-v7.3'); # added by Lu Cheng, 08.06.2012
#         else
#             tietue.proportionsIt = proportionsIt;
#             tietue.pvalue = uskottavuus; # Added by Jing
#             tietue.mixtureType = 'admix';
#             tietue.admixnpops = npops;
#     #         save([pathname filename], 'tietue');
#             save([pathname filename], 'tietue', '-v7.3'); # added by Lu Cheng, 08.06.2012
#         end
#     end
}


# function [npops] = poistaLiianPienet(npops, rowsFromInd, alaraja)
# % Muokkaa tulokset muotoon, jossa outlier yksilöt on
# % poistettu. Tarkalleen ottaen poistaa ne populaatiot, 
# % joissa on vähemmän kuin 'alaraja':n verran yksilöit?

# global PARTITION;
# global COUNTS;
# global SUMCOUNTS;

# popSize=zeros(1,npops);
# for i=1:npops
#     popSize(i)=length(find(PARTITION==i));
# end
# miniPops = find(popSize<alaraja);

# if length(miniPops)==0
#     return;
# end

# outliers = [];
# for pop = miniPops
#     inds = find(PARTITION==pop);
#     disp('Removed individuals: ');
#     disp(num2str(inds));
#     outliers = [outliers; inds];
# end

# ninds = length(PARTITION);
# PARTITION(outliers) = 0;
# korit = unique(PARTITION(find(PARTITION>0)));
# for n=1:length(korit)
#     kori = korit(n);
#     yksilot = find(PARTITION==kori);
#     PARTITION(yksilot) = n;
# end
# COUNTS(:,:,miniPops) = [];
# SUMCOUNTS(miniPops,:) = [];

# npops = npops-length(miniPops);

# %------------------------------------------------------------------------

# function clearGlobalVars

# global COUNTS; COUNTS = [];
# global SUMCOUNTS; SUMCOUNTS = [];
# global PARTITION; PARTITION = [];
# global POP_LOGML; POP_LOGML = [];

# %--------------------------------------------------------


# function allFreqs = computeAllFreqs2(noalle)
# % Lisää a priori jokaista alleelia
# % joka populaation joka lokukseen j 1/noalle(j) verran.

# global COUNTS;
# global SUMCOUNTS;

# max_noalle = size(COUNTS,1);
# nloci = size(COUNTS,2);
# npops = size(COUNTS,3);

# sumCounts = SUMCOUNTS+ones(size(SUMCOUNTS));
# sumCounts = reshape(sumCounts', [1, nloci, npops]);
# sumCounts = repmat(sumCounts, [max_noalle, 1 1]);

# prioriAlleelit = zeros(max_noalle,nloci);
# for j=1:nloci
#     prioriAlleelit(1:noalle(j),j) = 1/noalle(j);
# end
# prioriAlleelit = repmat(prioriAlleelit, [1,1,npops]);
# counts = COUNTS + prioriAlleelit;
# allFreqs = counts./sumCounts;


# function allfreqs = simulateAllFreqs(noalle)
# % Lisää jokaista alleelia joka populaation joka lokukseen j 1/noalle(j)
# % verran. Näin saatuja counts:eja vastaavista Dirichlet-jakaumista
# % simuloidaan arvot populaatioiden alleelifrekvensseille.

# global COUNTS;

# max_noalle = size(COUNTS,1);
# nloci = size(COUNTS,2);
# npops = size(COUNTS,3);

# prioriAlleelit = zeros(max_noalle,nloci);
# for j=1:nloci
#     prioriAlleelit(1:noalle(j),j) = 1/noalle(j);
# end
# prioriAlleelit = repmat(prioriAlleelit, [1,1,npops]);
# counts = COUNTS + prioriAlleelit;
# allfreqs = zeros(size(counts));

# for i=1:npops
#     for j=1:nloci
#         simuloidut = randdir(counts(1:noalle(j),j,i) , noalle(j));
#         allfreqs(1:noalle(j),j,i) = simuloidut;
#     end
# end

# %--------------------------------------------------------------------------


# function refData = simulateIndividuals(n,rowsFromInd,allfreqs,pop, missing_level)
# % simulate n individuals from population pop, such that approximately
# % proportion "missing_level" of the alleles are present. 

# nloci = size(allfreqs,2);

# refData = zeros(n*rowsFromInd,nloci);
# counter = 1;  % which row will be generated next.

# for ind = 1:n
#     for loc = 1:nloci
#         for k=0:rowsFromInd-1
#             if rand<missing_level
#                 refData(counter+k,loc) = simuloiAlleeli(allfreqs,pop,loc);
#             else
#                 refData(counter+k,loc) = -999;
#             end
#         end
#     end
#     counter = counter+rowsFromInd;
# end

# function all = simuloiAlleeli(allfreqs,pop,loc)
# % Simuloi populaation pop lokukseen loc alleelin.
# freqs = allfreqs(:,loc,pop);
# cumsumma = cumsum(freqs);
# arvo = rand;
# isommat = find(cumsumma>arvo);
# all = min(isommat);


# %--------------------------------------------------------------------------


# function omaFreqs = computePersonalAllFreqs(ind, data, allFreqs, rowsFromInd)
# % Laskee npops*(rowsFromInd*nloci) taulukon, jonka kutakin saraketta
# % vastaa yksilön ind alleeli. Eri rivit ovat alleelin alkuperäfrekvenssit
# % eri populaatioissa. Jos yksilölt?puuttuu jokin alleeli, niin vastaavaan
# % kohtaa tulee sarake ykkösi?

# global COUNTS;
# nloci = size(COUNTS,2);
# npops = size(COUNTS,3);

# rows = data(computeRows(rowsFromInd, ind, 1),:);

# omaFreqs = zeros(npops, (rowsFromInd*nloci));
# pointer = 1;
# for loc=1:size(rows,2)
#     for all=1:size(rows,1)
#         if rows(all,loc)>=0
#             try,
#             omaFreqs(:,pointer) = ...
#                 reshape(allFreqs(rows(all,loc),loc,:), [npops,1]);
#             catch
#                 a=0;
#             end
#         else
#             omaFreqs(:,pointer) = ones(npops,1);
#         end
#         pointer = pointer+1;
#     end
# end