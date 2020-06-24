#' @title Clustering of individuals
#' @param tietue Record
#' @param format Format of the data ("BAPS", "GenePop" or "Preprocessed")
#' @param savePreProcessed Save the pre-processed data?
#' @export
greedyMix <- function(
	tietue,
	format           = NULL,
	savePreProcessed = NULL,
	filePreProcessed = NULL
) {
	# ASK: graphical components. Remove?
	# check whether fixed k mode is selected
	# h0 <- findobj('Tag','fixk_menu')
	# fixedK = get(h0, 'userdata');

	# if fixedK
	# 	if ~(fixKWarning == 1) % call function fixKWarning
	# 		return
	# 	end
	# end

	# % check whether partition compare mode is selected
	# h1 = findobj('Tag','partitioncompare_menu');
	# partitionCompare = get(h1, 'userdata');

	if (is(tietue, "list") | is(tietue, "character")) {
		# ----------------------------------------------------------------------
		# Defining type of file
		# ----------------------------------------------------------------------
		if (is.null(format)) {
			input_type <- inputdlg(
				paste(
					'Specify the format of your data:\n',
					'1) BAPS-format\n',
					'2) GenePop-format\n',
					'3) Preprocessed data\n'
				)
			)
			# Converting from number into name
			input_type_name <- switch(
				input_type, 'BAPS-format', 'GenePop-format', 'Preprocessed data'
			)
		} else {
			input_type_name <- paste0(format, "-format")
		}
		# ----------------------------------------------------------------------
		# Treating BAPS-formatted files
		# ----------------------------------------------------------------------
		if (length(input_type_name) == 0) {
			stop('Invalid alternative')
		} else if (input_type_name == 'BAPS-format') {
			if (!is(tietue, "character")) {
				pathname_filename <- uigetfile("*.txt", "Loaddata in BAPS-format")
			} else {
				pathname_filename <- tietue
			}

			# ASK: remove?
			# if ~isempty(partitionCompare)
			# 	fprintf(1,'Data: %s\n',[pathname filename]);
			# end

			data <- read.delim(pathname_filename) # ASK: what is the delimiter?
			# ninds <- testaaOnkoKunnollinenBapsData(data)  #TESTAUS # TODO: trans
			# if (ninds == 0) stop('Incorrect Data-file')

			# ASK: remove?
			# h0 = findobj('Tag','filename1_text');
			# set(h0,'String',filename); clear h0;
			cat(
				'When using data which are in BAPS-format,',
				'you can specify the sampling populations of the',
				'individuals by giving two additional files:',
				'one containing the names of the populations,',
				'the other containing the indices of the first',
				'individuals of the populations.'
			)
			input_pops <- inputdlg(
				prompt = 'Do you wish to specify the sampling populations? [y/N]',
				definput = 'N'
			)
			if (tolower(input_pops) %in% c('yes', 'y')) {
				popfile    <- uigetfile('*.txt', 'Load population names')
				kysyToinen <- ifelse(popfile$name == 0, 0, 1)
				if (kysyToinen == 1) {
					indicesfile <- uigetfile('*.txt', 'Load population indices')
					if (indicesfile == 0) {
						popnames = ""
					} else {
						# popnames = initPopNames([namepath namefile],[indicespath indicesfile]) # TODO: translate this fun
					}
				} else {
					popnames <- ""
				}
			} else {
				popnames <- ""
			}

			# [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data); # TODO: translate this function
			# [Z,dist] = newGetDistances(data,rowsFromInd); # TODO: translate
			if (is.null(savePreProcessed)) {
				save_preproc <- questdlg(
					quest = 'Do you wish to save pre-processed data?',
					dlgtitle = 'Save pre-processed data?',
					defbtn = 'y'
				)
			} else {
				save_preproc <- savePreProcessed
			}
			if (save_preproc %in% c('y', 'yes', TRUE)) {
				file_out      <- uiputfile('.rda','Save pre-processed data as')
				kokonimi      <- paste0(file_out$path, file_out$name)
				c             <- list()
				c$data        <- data
				c$rowsFromInd <- rowsFromInd
				c$alleleCodes <- alleleCodes
				c$noalle      <- noalle
				c$adjprior    <- adjprior
				c$priorTerm   <- priorTerm
				c$dist        <- dist
				c$popnames    <- popnames
				c$Z           <- Z
				save(c, file = kokonimi)
				rm(c)
			}
		# ----------------------------------------------------------------------
		# Treating GenePop-formatted files
		# ----------------------------------------------------------------------
		} else if (input_type_name == 'GenePop-format') {
			if (!is(tietue, "character")) {
				filename_pathname <- uigetfile(
					filter = '*.txt',
					title = 'Load data in GenePop-format'
				)
				if (filename_pathname$name == 0) stop("No name provided")
			} else {
				filename_pathname <- tietue
			}

			# ASK: remove?
			# if ~isempty(partitionCompare)
			# 	fprintf(1,'Data: %s\n',[pathname filename]);
			# end

			# kunnossa = testaaGenePopData([pathname filename]); # TODO: trans
			# if (kunnossa == 0) stop("testaaGenePopData returned 0")
			# [data,popnames]=lueGenePopData([pathname filename]); # TODO: trans

# 			h0 = findobj('Tag','filename1_text');
# 			set(h0,'String',filename); clear h0;

# 			[data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data); # TODO:trans
# 			[Z,dist] = newGetDistances(data,rowsFromInd); # TODO: trans
			if (is.null(savePreProcessed)) {
				save_preproc <- questdlg(
					quest = 'Do you wish to save pre-processed data?',
					dlgtitle = 'Save pre-processed data?',
					defbtn = 'y'
				)
			} else {
				save_preproc <- savePreProcessed
			}
			if (save_preproc %in% c('y', 'Yes', TRUE)) {
				file_out      <- uiputfile('.rda','Save pre-processed data as')
				kokonimi      <- paste0(file_out$path, file_out$name)
				# FIXME: translate functions above so the objects below exist
				c$data        <- data
				c$rowsFromInd <- rowsFromInd
				c$alleleCodes <- alleleCodes
				c$noalle      <- noalle
				c$adjprior    <- adjprior
				c$priorTerm   <- priorTerm
				c$dist        <- dist
				c$popnames    <- popnames
				c$Z           <- Z
				save(c, file = kokonimi)
				rm(c)
			}
		# ----------------------------------------------------------------------
		# Handling Pre-processed data
		# ----------------------------------------------------------------------
		} else if (input_type_name == 'Preprocessed data') {
			file_in <- uigetfile(
				filter = '*.txt',
				title = 'Load pre-processed data in GenePop-format'
			)
			if (file_in$name == 0) stop("No name provided")

			# ASK: remove?
			# h0 = findobj('Tag','filename1_text');
			# set(h0,'String',filename); clear h0;
			# if ~isempty(partitionCompare)
			# 	fprintf(1,'Data: %s\n',[pathname filename]);
			# end

			struct_array <- readRDS(paste0(file_in$path, file_in$name))
			if (isfield(struct_array,'c')) { # Matlab versio
				c <- struct_array$c
				if (!isfield(c,'dist')) stop('Incorrect file format')
			} else if (isfield(struct_array,'dist')) { #Mideva versio
				c <- struct_array
			} else {
				stop('Incorrect file format')
			}
			data        <- double(c$data)
			rowsFromInd <- c$rowsFromInd
			alleleCodes <- c$alleleCodes
			noalle      <- c$noalle
			adjprior    <- c$adjprior
			priorTerm   <- c$priorTerm
			dist        <- c$dist
			popnames    <- c$popnames
			Z           <- c$Z
			rm(c)
		}
	} else {
		data        <- double(tietue$data)
		rowsFromInd <- tietue$rowsFromInd
		alleleCodes <- tietue$alleleCodes
		noalle      <- tietue$noalle
		adjprior    <- tietue$adjprior
		priorTerm   <- tietue$priorTerm
		dist        <- tietue$dist
		popnames    <- tietue$popnames
		Z           <- tietue$Z
		rm(tietue)
	}

	# ==========================================================================
	# Declaring global variables
	# ==========================================================================
	PARTITION       <- vector()
	COUNTS          <- vector()
	SUMCOUNTS       <- vector()
	POP_LOGML       <- vector()
	clearGlobalVars <- vector()
	# ==========================================================================

	c$data        <- data
	c$noalle      <- noalle
	c$adjprior    <- adjprior
	c$priorTerm   <- priorTerm
	c$dist        <- dist
	c$Z           <- Z
	c$rowsFromInd <- rowsFromInd

	ninds  <- length(unique(data[, ncol(data)]))
	ekat   <- t(seq(1, ninds, rowsFromInd) * rowsFromInd)
	c$rows <- c(ekat, ekat + rowsFromInd - 1)

	# partition compare
	if (!is.null(partitionCompare)) {
		nsamplingunits <- size(c$rows, 1)
		partitions     <- partitionCompare$partitions
		npartitions    <- size(partitions, 2)
		partitionLogml <- zeros(1, npartitions)
		for (i in seq_len(npartitions)) {
			# number of unique partition lables
			npops <- length(unique(partitions[, i]))

			partitionInd <- zeros(ninds * rowsFromInd, 1)
			partitionSample <- partitions[, i]
			for (j in seq_len(nsamplingunits)) {
				partitionInd[c$rows[j, 1]:c$rows[j, 2]] <- partitionSample[j]
			}
			# partitionLogml[i] = initialCounts(
			# 	partitionInd,
			# 	data[, seq_len(end - 1)],
			# 	npops,
			# 	c$rows,
			# 	noalle,
			# 	adjprior
			# ) #TODO translate
		}
		# return the logml result
		partitionCompare$logmls <- partitionLogml
		# set(h1, 'userdata', partitionCompare) # ASK remove?
		return()
	}

	# ASK remove (graphical part)?
	# if (fixedK) {
	# 	#logml_npops_partitionSummary <- indMix_fixK(c) # TODO translate?
	# } else {
	# 	#logml_npops_partitionSummary <- indMix(c) # TODO translate?
	# }
	# if (logml_npops_partitionSummary$logml == 1) return()

	data <- data[, seq_len(ncol(data) - 1)]

	# ASK: remove?
	# h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
	# h0 = findobj('Tag','filename2_text');
	# outp = get(h0,'String');

	# changesInLogml <- writeMixtureInfo(
	# 	logml, rowsFromInd, data, adjprior, priorTerm, outp, inp,
	# 	popnames, fixedK
	# ) # TODO translate

	# viewMixPartition(PARTITION, popnames) # TODO translate function

	talle <- questdlg(
		quest = paste(
			'Do you want to save the mixture populations',
			'so that you can use them later in admixture analysis?'
		),
		dlgtitle = 'Save results?',
		defbtn = 'y'
	)
	if (talle %in% c('Yes', 'y')) {
		filename_pathname <- uiputfile('.mat','Save results as')

	# ==========================================================================
		cond <- (sum(filename_pathname$name) == 0) |
		(sum(filename_pathname$path) == 0)
		if (cond) {
			# Cancel was pressed
			return()
		} else {
			# copy 'baps4_output.baps' into the text file with the same name.
			if (file.exists('baps4_output.baps')) {
				file.copy(
					from = 'baps4_output.baps',
					to = paste0(
						filename_pathname$path, filename_pathname$name, '.txt'
					)
				)
				file.remove('baps4_output.baps')
			}
		}
	# ==========================================================================

		c$PARTITION      <- PARTITION
		c$COUNTS         <- COUNTS
		c$SUMCOUNTS      <- SUMCOUNTS
		c$alleleCodes    <- alleleCodes
		c$adjprior       <- adjprior
		c$popnames       <- popnames
		c$rowsFromInd    <- rowsFromInd
		c$data           <- data
		c$npops          <- npops
		c$noalle         <- noalle
		c$mixtureType    <- 'mix'
		c$logml          <- logml
		c$changesInLogml <- changesInLogml

		save(c, file = paste0(filename_pathname$path, filename_pathname$name))
	} else {
		if (file.exists('baps4_output.baps')) file.remove('baps4_output.baps')
	}
}

# 	%-------------------------------------------------------------------------------------
# 	%-------------------------------------------------------------------------------------

# 	function clearGlobalVars

# 	global COUNTS; COUNTS = [];
# 	global SUMCOUNTS; SUMCOUNTS = [];
# 	global PARTITION; PARTITION = [];
# 	global POP_LOGML; POP_LOGML = [];


# 	%-------------------------------------------------------------------------------------


# 	function rows = computeRows(rowsFromInd, inds, ninds)
# 	% On annettu yksil�t inds. Funktio palauttaa vektorin, joka
# 	% sis�lt�� niiden rivien numerot, jotka sis�lt�v�t yksil�iden
# 	% dataa.

# 	rows = inds(:, ones(1,rowsFromInd));
# 	rows = rows*rowsFromInd;
# 	miinus = repmat(rowsFromInd-1 : -1 : 0, [ninds 1]);
# 	rows = rows - miinus;
# 	rows = reshape(rows', [1,rowsFromInd*ninds]);


# 	%--------------------------------------------------------------------------


# 	function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
# 	% Tiedet��n, ett?annettu logml on isompi kuin huonoin arvo
# 	% partitionSummary taulukossa. Jos partitionSummary:ss?ei viel?ole
# 	% annettua logml arvoa, niin lis�t��n worstIndex:in kohtaan uusi logml ja
# 	% nykyist?partitiota vastaava nclusters:in arvo. Muutoin ei tehd?mit��n.

# 	apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
# 	if isempty(apu)
# 		% Nyt l�ydetty partitio ei ole viel?kirjattuna summaryyn.
# 		global PARTITION;
# 		npops = length(unique(PARTITION));
# 		partitionSummary(worstIndex,1) = npops;
# 		partitionSummary(worstIndex,2) = logml;
# 		added = 1;
# 	else
# 		added = 0;
# 	end


# 	%--------------------------------------------------------------------------


# 	function [suurin, i2] = arvoSeuraavaTila(muutokset, logml)
# 	% Suorittaa yksil�n seuraavan tilan arvonnan

# 	y = logml + muutokset;  % siirron j�lkeiset logml:t
# 	y = y - max(y);
# 	y = exp(y);
# 	summa = sum(y);
# 	y = y/summa;
# 	y = cumsum(y);

# 	i2 = rand_disc(y);   % uusi kori
# 	suurin = muutokset(i2);


# 	%--------------------------------------------------------------------------------------


# 	function svar=rand_disc(CDF)
# 	%returns an index of a value from a discrete distribution using inversion method
# 	slump=rand;
# 	har=find(CDF>slump);
# 	svar=har(1);


# 	%-------------------------------------------------------------------------------------


# 	function updateGlobalVariables(ind, i2, rowsFromInd, diffInCounts, ...
# 		adjprior, priorTerm)
# 	% Suorittaa globaalien muuttujien muutokset, kun yksil?ind
# 	% on siirret��n koriin i2.

# 	global PARTITION;
# 	global COUNTS;
# 	global SUMCOUNTS;
# 	global POP_LOGML;

# 	i1 = PARTITION(ind);
# 	PARTITION(ind)=i2;

# 	COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
# 	COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

# 	POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


# 	%---------------------------------------------------------------------------------


# 	function updateGlobalVariables2( ...
# 		i1, i2, rowsFromInd, diffInCounts, adjprior, priorTerm);
# 	% Suorittaa globaalien muuttujien muutokset, kun kaikki
# 	% korissa i1 olevat yksil�t siirret��n koriin i2.

# 	global PARTITION;
# 	global COUNTS;
# 	global SUMCOUNTS;
# 	global POP_LOGML;

# 	inds = find(PARTITION==i1);
# 	PARTITION(inds) = i2;

# 	COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
# 	COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

# 	POP_LOGML(i1) = 0;
# 	POP_LOGML(i2) = computePopulationLogml(i2, adjprior, priorTerm);


# 	%------------------------------------------------------------------------------------


# 	function updateGlobalVariables3(muuttuvat, rowsFromInd, diffInCounts, ...
# 		adjprior, priorTerm, i2);
# 	% Suorittaa globaalien muuttujien p�ivitykset, kun yksil�t 'muuttuvat'
# 	% siirret��n koriin i2. Ennen siirtoa yksil�iden on kuuluttava samaan
# 	% koriin.

# 	global PARTITION;
# 	global COUNTS;
# 	global SUMCOUNTS;
# 	global POP_LOGML;

# 	i1 = PARTITION(muuttuvat(1));
# 	PARTITION(muuttuvat) = i2;

# 	COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
# 	COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

# 	POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);


# 	%----------------------------------------------------------------------


# 	function inds = returnInOrder(inds, pop, rowsFromInd, data, adjprior, priorTerm)
# 	% Palauttaa yksil�t j�rjestyksess?siten, ett?ensimm�isen?on
# 	% se, jonka poistaminen populaatiosta pop nostaisi logml:n
# 	% arvoa eniten.

# 	global COUNTS;      global SUMCOUNTS;
# 	ninds = length(inds);
# 	apuTaulu = [inds, zeros(ninds,1)];

# 	for i=1:ninds
# 		ind = inds(i);
# 		rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
# 		diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
# 		diffInSumCounts = sum(diffInCounts);

# 		COUNTS(:,:,pop) = COUNTS(:,:,pop)-diffInCounts;
# 		SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)-diffInSumCounts;
# 		apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
# 		COUNTS(:,:,pop) = COUNTS(:,:,pop)+diffInCounts;
# 		SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)+diffInSumCounts;
# 	end
# 	apuTaulu = sortrows(apuTaulu,2);
# 	inds = apuTaulu(ninds:-1:1,1);

# 	%------------------------------------------------------------------------------------


# 	function [muutokset, diffInCounts] = ...
# 		laskeMuutokset(ind, rowsFromInd, data, adjprior, priorTerm)
# 	% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
# 	% muutos logml:ss? mik�li yksil?ind siirret��n koriin i.
# 	% diffInCounts on poistettava COUNTS:in siivusta i1 ja lis�tt�v?
# 	% COUNTS:in siivuun i2, mik�li muutos toteutetaan.

# 	global COUNTS;      global SUMCOUNTS;
# 	global PARTITION;   global POP_LOGML;
# 	npops = size(COUNTS,3);
# 	muutokset = zeros(npops,1);

# 	i1 = PARTITION(ind);
# 	i1_logml = POP_LOGML(i1);

# 	rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
# 	diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
# 	diffInSumCounts = sum(diffInCounts);

# 	COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
# 	new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
# 	COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

# 	i2 = [1:i1-1 , i1+1:npops];
# 	i2_logml = POP_LOGML(i2);

# 	COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
# 	new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
# 	COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

# 	muutokset(i2) = new_i1_logml - i1_logml ...
# 		+ new_i2_logml - i2_logml;


# 	%------------------------------------------------------------------------------------


# 	function [muutokset, diffInCounts] = laskeMuutokset2( ...
# 		i1, rowsFromInd, data, adjprior, priorTerm);
# 	% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik?olisi
# 	% muutos logml:ss? mik�li korin i1 kaikki yksil�t siirret��n
# 	% koriin i.

# 	global COUNTS;      global SUMCOUNTS;
# 	global PARTITION;   global POP_LOGML;
# 	npops = size(COUNTS,3);
# 	muutokset = zeros(npops,1);

# 	i1_logml = POP_LOGML(i1);

# 	inds = find(PARTITION==i1);
# 	ninds = length(inds);

# 	if ninds==0
# 		diffInCounts = zeros(size(COUNTS,1), size(COUNTS,2));
# 		return;
# 	end

# 	rows = computeRows(rowsFromInd, inds, ninds);

# 	diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
# 	diffInSumCounts = sum(diffInCounts);

# 	COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
# 	new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
# 	COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
# 	SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

# 	i2 = [1:i1-1 , i1+1:npops];
# 	i2_logml = POP_LOGML(i2);

# 	COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
# 	new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
# 	COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
# 	SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

# 	muutokset(i2) = new_i1_logml - i1_logml ...
# 		+ new_i2_logml - i2_logml;



# 	%------------------------------------------------------------------------------------


# 	function muutokset = laskeMuutokset3(T2, inds2, rowsFromInd, ...
# 		data, adjprior, priorTerm, i1)
# 	% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
# 	% kertoo, mik?olisi muutos logml:ss? jos populaation i1 osapopulaatio
# 	% inds2(find(T2==i)) siirret��n koriin j.

# 	global COUNTS;      global SUMCOUNTS;
# 	global PARTITION;   global POP_LOGML;
# 	npops = size(COUNTS,3);
# 	npops2 = length(unique(T2));
# 	muutokset = zeros(npops2, npops);

# 	i1_logml = POP_LOGML(i1);

# 	for pop2 = 1:npops2
# 		inds = inds2(find(T2==pop2));
# 		ninds = length(inds);
# 		if ninds>0
# 			rows = computeRows(rowsFromInd, inds, ninds);
# 			diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
# 			diffInSumCounts = sum(diffInCounts);

# 			COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
# 			SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
# 			new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
# 			COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
# 			SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

# 			i2 = [1:i1-1 , i1+1:npops];
# 			i2_logml = POP_LOGML(i2)';

# 			COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
# 			SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
# 			new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
# 			COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
# 			SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

# 			muutokset(pop2,i2) = new_i1_logml - i1_logml ...
# 				+ new_i2_logml - i2_logml;
# 		end
# 	end


# 	%------------------------------------------------------------------------------------

# 	function muutokset = laskeMuutokset5(inds, rowsFromInd, data, adjprior, ...
# 		priorTerm, i1, i2)

# 	% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik?olisi
# 	% muutos logml:ss? mik�li yksil?i vaihtaisi koria i1:n ja i2:n v�lill?

# 	global COUNTS;      global SUMCOUNTS;
# 	global PARTITION;   global POP_LOGML;

# 	ninds = length(inds);
# 	muutokset = zeros(ninds,1);

# 	i1_logml = POP_LOGML(i1);
# 	i2_logml = POP_LOGML(i2);

# 	for i = 1:ninds
# 		ind = inds(i);
# 		if PARTITION(ind)==i1
# 			pop1 = i1;  %mist?
# 			pop2 = i2;  %mihin
# 		else
# 			pop1 = i2;
# 			pop2 = i1;
# 		end
# 		rows = (ind-1)*rowsFromInd+1 : ind*rowsFromInd;
# 		diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
# 		diffInSumCounts = sum(diffInCounts);

# 		COUNTS(:,:,pop1) = COUNTS(:,:,pop1)-diffInCounts;
# 		SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)-diffInSumCounts;
# 		COUNTS(:,:,pop2) = COUNTS(:,:,pop2)+diffInCounts;
# 		SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)+diffInSumCounts;

# 		PARTITION(ind) = pop2;

# 		new_logmls = computePopulationLogml([i1 i2], adjprior, priorTerm);

# 		muutokset(i) = sum(new_logmls);

# 		COUNTS(:,:,pop1) = COUNTS(:,:,pop1)+diffInCounts;
# 		SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)+diffInSumCounts;
# 		COUNTS(:,:,pop2) = COUNTS(:,:,pop2)-diffInCounts;
# 		SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)-diffInSumCounts;

# 		PARTITION(ind) = pop1;
# 	end

# 	muutokset = muutokset - i1_logml - i2_logml;

# 	%--------------------------------------------------------------------------



# 	function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
# 	% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
# 	% lukum��r�t (vastaavasti kuin COUNTS:issa), jotka ovat data:n
# 	% riveill?rows.

# 	diffInCounts = zeros(max_noalle, nloci);
# 	for i=rows
# 		row = data(i,:);
# 		notEmpty = find(row>=0);

# 		if length(notEmpty)>0
# 			diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
# 				diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
# 		end
# 	end



# 	%------------------------------------------------------------------------------------


# 	function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
# 	% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
# 	% logml:t koreille, jotka on m��ritelty pops-muuttujalla.

# 	global COUNTS;
# 	global SUMCOUNTS;
# 	x = size(COUNTS,1);
# 	y = size(COUNTS,2);
# 	z = length(pops);

# 	popLogml = ...
# 		squeeze(sum(sum(reshape(...
# 		gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
# 		,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;


# 	%-----------------------------------------------------------------------------------


# 	function npops = poistaTyhjatPopulaatiot(npops)
# 	% Poistaa tyhjentyneet populaatiot COUNTS:ista ja
# 	% SUMCOUNTS:ista. P�ivitt�� npops:in ja PARTITION:in.

# 	global COUNTS;
# 	global SUMCOUNTS;
# 	global PARTITION;

# 	notEmpty = find(any(SUMCOUNTS,2));
# 	COUNTS = COUNTS(:,:,notEmpty);
# 	SUMCOUNTS = SUMCOUNTS(notEmpty,:);

# 	for n=1:length(notEmpty)
# 		apu = find(PARTITION==notEmpty(n));
# 		PARTITION(apu)=n;
# 	end
# 	npops = length(notEmpty);


# 	%----------------------------------------------------------------------------------
# 	%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

# 	function initial_partition=admixture_initialization(data_matrix,nclusters,Z)
# 	size_data=size(data_matrix);
# 	nloci=size_data(2)-1;
# 	n=max(data_matrix(:,end));
# 	T=cluster_own(Z,nclusters);
# 	initial_partition=zeros(size_data(1),1);
# 	for i=1:n
# 		kori=T(i);
# 		here=find(data_matrix(:,end)==i);
# 		for j=1:length(here)
# 			initial_partition(here(j),1)=kori;
# 		end
# 	end

# 	function T = cluster_own(Z,nclust)
# 	true=logical(1);
# 	false=logical(0);
# 	maxclust = nclust;
# 	% Start of algorithm
# 	m = size(Z,1)+1;
# 	T = zeros(m,1);
# 	   % maximum number of clusters based on inconsistency
# 	   if m <= maxclust
# 		  T = (1:m)';
# 	   elseif maxclust==1
# 		  T = ones(m,1);
# 	   else
# 		  clsnum = 1;
# 		  for k = (m-maxclust+1):(m-1)
# 			 i = Z(k,1); % left tree
# 			 if i <= m % original node, no leafs
# 				T(i) = clsnum;
# 				clsnum = clsnum + 1;
# 			 elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
# 				T = clusternum(Z, T, i-m, clsnum);
# 				clsnum = clsnum + 1;
# 			 end
# 			 i = Z(k,2); % right tree
# 			 if i <= m  % original node, no leafs
# 				T(i) = clsnum;
# 				clsnum = clsnum + 1;
# 			 elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
# 				T = clusternum(Z, T, i-m, clsnum);
# 				clsnum = clsnum + 1;
# 			 end
# 		  end
# 	   end

# 	function T = clusternum(X, T, k, c)
# 	m = size(X,1)+1;
# 	while(~isempty(k))
# 	   % Get the children of nodes at this level
# 	   children = X(k,1:2);
# 	   children = children(:);

# 	   % Assign this node number to leaf children
# 	   t = (children<=m);
# 	   T(children(t)) = c;

# 	   % Move to next level
# 	   k = children(~t) - m;
# 	end


# 	%---------------------------------------------------------------------------------------


# 	function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
# 		handleData(raw_data)
# 	% Alkuper�isen datan viimeinen sarake kertoo, milt?yksil�lt?
# 	% kyseinen rivi on per�isin. Funktio tutkii ensin, ett?montako
# 	% rivi?maksimissaan on per�isin yhdelt?yksil�lt? jolloin saadaan
# 	% tiet�� onko kyseess?haploidi, diploidi jne... T�m�n j�lkeen funktio
# 	% lis�� tyhji?rivej?niille yksil�ille, joilta on per�isin v�hemm�n
# 	% rivej?kuin maksimim��r?
# 	%   Mik�li jonkin alleelin koodi on =0, funktio muuttaa t�m�n alleelin
# 	% koodi pienimm�ksi koodiksi, joka isompi kuin mik��n k�yt�ss?oleva koodi.
# 	% T�m�n j�lkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
# 	% koodit saavat arvoja v�lill?1,...,noalle(j).
# 	data = raw_data;
# 	nloci=size(raw_data,2)-1;

# 	dataApu = data(:,1:nloci);
# 	nollat = find(dataApu==0);
# 	if ~isempty(nollat)
# 	   isoinAlleeli = max(max(dataApu));
# 	   dataApu(nollat) = isoinAlleeli+1;
# 	   data(:,1:nloci) = dataApu;
# 	end
# 	dataApu = []; nollat = []; isoinAlleeli = [];

# 	noalle=zeros(1,nloci);
# 	alleelitLokuksessa = cell(nloci,1);
# 	for i=1:nloci
# 		alleelitLokuksessaI = unique(data(:,i));
# 		alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
# 		noalle(i) = length(alleelitLokuksessa{i,1});
# 	end
# 	alleleCodes = zeros(max(noalle),nloci);
# 	for i=1:nloci
# 		alleelitLokuksessaI = alleelitLokuksessa{i,1};
# 		puuttuvia = max(noalle)-length(alleelitLokuksessaI);
# 		alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
# 	end

# 	for loc = 1:nloci
# 		for all = 1:noalle(loc)
# 			data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
# 		end;
# 	end;

# 	nind = max(data(:,end));
# 	nrows = size(data,1);
# 	ncols = size(data,2);
# 	rowsFromInd = zeros(nind,1);
# 	for i=1:nind
# 		rowsFromInd(i) = length(find(data(:,end)==i));
# 	end
# 	maxRowsFromInd = max(rowsFromInd);
# 	a = -999;
# 	emptyRow = repmat(a, 1, ncols);
# 	lessThanMax = find(rowsFromInd < maxRowsFromInd);
# 	missingRows = maxRowsFromInd*nind - nrows;
# 	data = [data; zeros(missingRows, ncols)];
# 	pointer = 1;
# 	for ind=lessThanMax'    %K�y l�pi ne yksil�t, joilta puuttuu rivej?
# 		miss = maxRowsFromInd-rowsFromInd(ind);  % T�lt?yksil�lt?puuttuvien lkm.
# 		for j=1:miss
# 			rowToBeAdded = emptyRow;
# 			rowToBeAdded(end) = ind;
# 			data(nrows+pointer, :) = rowToBeAdded;
# 			pointer = pointer+1;
# 		end
# 	end
# 	data = sortrows(data, ncols);   % Sorttaa yksil�iden mukaisesti
# 	newData = data;
# 	rowsFromInd = maxRowsFromInd;

# 	adjprior = zeros(max(noalle),nloci);
# 	priorTerm = 0;
# 	for j=1:nloci
# 		adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
# 		priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
# 	end


# 	%----------------------------------------------------------------------------------------


# 	function [Z, dist] = newGetDistances(data, rowsFromInd)

# 	ninds = max(data(:,end));
# 	nloci = size(data,2)-1;
# 	riviLkm = nchoosek(ninds,2);

# 	empties = find(data<0);
# 	data(empties)=0;
# 	data = uint8(data);   % max(noalle) oltava <256

# 	pariTaulu = zeros(riviLkm,2);
# 	aPointer=1;
# 	for a=1:ninds-1
# 		pariTaulu(aPointer:aPointer+ninds-1-a,1) = ones(ninds-a,1)*a;
# 		pariTaulu(aPointer:aPointer+ninds-1-a,2) = (a+1:ninds)';
# 		aPointer = aPointer+ninds-a;
# 	end

# 	eka = pariTaulu(:,ones(1,rowsFromInd));
# 	eka = eka * rowsFromInd;
# 	miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
# 	eka = eka - miinus;

# 	toka = pariTaulu(:,ones(1,rowsFromInd)*2);
# 	toka = toka * rowsFromInd;
# 	toka = toka - miinus;

# 	%eka = uint16(eka);
# 	%toka = uint16(toka);

# 	summa = zeros(riviLkm,1);
# 	vertailuja = zeros(riviLkm,1);

# 	clear pariTaulu; clear miinus;

# 	x = zeros(size(eka));    x = uint8(x);
# 	y = zeros(size(toka));   y = uint8(y);

# 	for j=1:nloci;

# 		for k=1:rowsFromInd
# 			x(:,k) = data(eka(:,k),j);
# 			y(:,k) = data(toka(:,k),j);
# 		end

# 		for a=1:rowsFromInd
# 			for b=1:rowsFromInd
# 				vertailutNyt = double(x(:,a)>0 & y(:,b)>0);
# 				vertailuja = vertailuja + vertailutNyt;
# 				lisays = (x(:,a)~=y(:,b) & vertailutNyt);
# 				summa = summa+double(lisays);
# 			end
# 		end
# 	end

# 	clear x;    clear y;   clear vertailutNyt;
# 	nollat = find(vertailuja==0);
# 	dist = zeros(length(vertailuja),1);
# 	dist(nollat) = 1;
# 	muut = find(vertailuja>0);
# 	dist(muut) = summa(muut)./vertailuja(muut);
# 	clear summa; clear vertailuja;

# 	Z = linkage(dist');


# 	%----------------------------------------------------------------------------------------


# 	function [Z, distances]=getDistances(data_matrix,nclusters)

# 	%finds initial admixture clustering solution with nclusters clusters, uses simple mean Hamming distance
# 	%gives partition in 8-bit format
# 	%allocates all alleles of a single individual into the same basket
# 	%data_matrix contains #Loci+1 columns, last column indicate whose alleles are placed in each row,
# 	%i.e. ranges from 1 to #individuals. For diploids there are 2 rows per individual, for haploids only a single row
# 	%missing values are indicated by zeros in the partition and by negative integers in the data_matrix.

# 	size_data=size(data_matrix);
# 	nloci=size_data(2)-1;
# 	n=max(data_matrix(:,end));
# 	distances=zeros(nchoosek(n,2),1);
# 	pointer=1;
# 	for i=1:n-1
# 		i_data=data_matrix(find(data_matrix(:,end)==i),1:nloci);
# 		for j=i+1:n
# 			d_ij=0;
# 			j_data=data_matrix(find(data_matrix(:,end)==j),1:nloci);
# 			vertailuja = 0;
# 			for k=1:size(i_data,1)
# 				for l=1:size(j_data,1)
# 					here_i=find(i_data(k,:)>=0);
# 					here_j=find(j_data(l,:)>=0);
# 					here_joint=intersect(here_i,here_j);
# 					vertailuja = vertailuja + length(here_joint);
# 					d_ij = d_ij + length(find(i_data(k,here_joint)~=j_data(l,here_joint)));
# 				end
# 			end
# 			d_ij = d_ij / vertailuja;
# 			distances(pointer)=d_ij;
# 			pointer=pointer+1;
# 		end
# 	end

# 	Z=linkage(distances');



# 	%----------------------------------------------------------------------------------------


# 	function Z = linkage(Y, method)
# 	[k, n] = size(Y);
# 	m = (1+sqrt(1+8*n))/2;
# 	if k ~= 1 | m ~= fix(m)
# 	  error('The first input has to match the output of the PDIST function in size.');
# 	end
# 	if nargin == 1 % set default switch to be 'co'
# 	   method = 'co';
# 	end
# 	method = lower(method(1:2)); % simplify the switch string.
# 	monotonic = 1;
# 	Z = zeros(m-1,3); % allocate the output matrix.
# 	N = zeros(1,2*m-1);
# 	N(1:m) = 1;
# 	n = m; % since m is changing, we need to save m in n.
# 	R = 1:n;
# 	for s = 1:(n-1)
# 	   X = Y;
# 	   [v, k] = min(X);
# 	   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
# 	   j = k - (i-1)*(m-i/2)+i;
# 	   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
# 	   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
# 	   U = [I1 I2 I3];
# 	   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
# 	   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

# 	   switch method
# 	   case 'si' %single linkage
# 		  Y(I) = min(Y(I),Y(J));
# 	   case 'av' % average linkage
# 		  Y(I) = Y(I) + Y(J);
# 	   case 'co' %complete linkage
# 		  Y(I) = max(Y(I),Y(J));
# 	   case 'ce' % centroid linkage
# 		  K = N(R(i))+N(R(j));
# 		  Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
# 	   case 'wa'
# 		  Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
# 		  N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
# 	   end
# 	   J = [J i*(m-(i+1)/2)-m+j];
# 	   Y(J) = []; % no need for the cluster information about j.

# 	   % update m, N, R
# 	   m = m-1;
# 	   N(n+s) = N(R(i)) + N(R(j));
# 	   R(i) = n+s;
# 	   R(j:(n-1))=R((j+1):n);
# 	end


# 	%-----------------------------------------------------------------------------------


# 	function popnames = initPopNames(nameFile, indexFile)
# 	%Palauttaa tyhj�n, mik�li nimitiedosto ja indeksitiedosto
# 	% eiv�t olleet yht?pitki?

# 	popnames = [];
# 	indices = load(indexFile);

# 	fid = fopen(nameFile);
# 	if fid == -1
# 		%File didn't exist
# 		msgbox('Loading of the population names was unsuccessful', ...
# 			'Error', 'error');
# 		return;
# 	end;
# 	line = fgetl(fid);
# 	counter = 1;
# 	while (line ~= -1) & ~isempty(line)
# 		names{counter} = line;
# 		line = fgetl(fid);
# 		counter = counter + 1;
# 	end;
# 	fclose(fid);

# 	if length(names) ~= length(indices)
# 		disp('The number of population names must be equal to the number of ');
# 		disp('entries in the file specifying indices of the first individuals of ');
# 		disp('each population.');
# 		return;
# 	end

# 	popnames = cell(length(names), 2);
# 	for i = 1:length(names)
# 		popnames{i,1} = names(i);
# 		popnames{i,2} = indices(i);
# 	end


# 	%-----------------------------------------------------------------------------------


# 	%--------------------------------------------------------------------------

# 	function logml = ...
# 		initialCounts(partition, data, npops, rows, noalle, adjprior)

# 	nloci=size(data,2);
# 	ninds = size(rows, 1);

# 	koot = rows(:,1) - rows(:,2) + 1;
# 	maxSize = max(koot);

# 	counts = zeros(max(noalle),nloci,npops);
# 	sumcounts = zeros(npops,nloci);
# 	for i=1:npops
# 		for j=1:nloci
# 			havainnotLokuksessa = find(partition==i & data(:,j)>=0);
# 			sumcounts(i,j) = length(havainnotLokuksessa);
# 			for k=1:noalle(j)
# 				alleleCode = k;
# 				N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
# 				counts(k,j,i) = N_ijk;
# 			end
# 		end
# 	end

# 	%initializeGammaln(ninds, maxSize, max(noalle));

# 	logml = laskeLoggis(counts, sumcounts, adjprior);


# 	%-----------------------------------------------------------------------


# 	function logml=computeLogml(counts, sumcounts, noalle, data, rowsFromInd)
# 	nloci = size(counts,2);
# 	npops = size(counts,3);
# 	adjnoalle = zeros(max(noalle),nloci);
# 	for j=1:nloci
# 		adjnoalle(1:noalle(j),j)=noalle(j);
# 		if (noalle(j)<max(noalle))
# 			adjnoalle(noalle(j)+1:end,j)=1;
# 		end
# 	end

# 	%logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
# 	%    - npops*sum(sum(gammaln(adjprior))) - ...
# 	%    sum(sum(gammaln(1+sumcounts)));
# 	%logml = logml2;

# 	global GAMMA_LN;
# 	rowsInG = size(data,1)+rowsFromInd;

# 	logml = sum(sum(sum(GAMMA_LN(counts+1 + repmat(rowsInG*(adjnoalle-1),[1 1 npops]))))) ...
# 		- npops*sum(sum(GAMMA_LN(1, adjnoalle))) ...
# 		-sum(sum(GAMMA_LN(sumcounts+1,1)));


# 	%--------------------------------------------------------------------------


# 	function initializeGammaln(ninds, rowsFromInd, maxAlleles)
# 	%Alustaa GAMMALN muuttujan s.e. GAMMALN(i,j)=gammaln((i-1) + 1/j)
# 	global GAMMA_LN;
# 	GAMMA_LN = zeros((1+ninds)*rowsFromInd, maxAlleles);
# 	for i=1:(ninds+1)*rowsFromInd
# 		for j=1:maxAlleles
# 			GAMMA_LN(i,j)=gammaln((i-1) + 1/j);
# 		end
# 	end

# 	%-------------------------------------------------------------------


# 	function changesInLogml = writeMixtureInfo(logml, rowsFromInd, data, adjprior, ...
# 		priorTerm, outPutFile, inputFile, partitionSummary, popnames, fixedK)

# 	global PARTITION;
# 	global COUNTS;
# 	global SUMCOUNTS;
# 	global LOGDIFF;
# 	changesInLogml = [];
# 	ninds = size(data,1)/rowsFromInd;
# 	npops =  size(COUNTS,3);
# 	names = (size(popnames,1) == ninds);    %Tarkistetaan ett?nimet viittaavat yksil�ihin

# 	if length(outPutFile)>0
# 		fid = fopen(outPutFile,'a');
# 	else
# 		fid = -1;
# 		diary('baps4_output.baps'); % save in text anyway.
# 	end

# 	dispLine;
# 	disp('RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:');
# 	disp(['Data file: ' inputFile]);
# 	disp(['Model: independent']);
# 	disp(['Number of clustered individuals: ' ownNum2Str(ninds)]);
# 	disp(['Number of groups in optimal partition: ' ownNum2Str(npops)]);
# 	disp(['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]);
# 	disp(' ');
# 	if (fid ~= -1)
# 		fprintf(fid,'%s \n', ['RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:']); fprintf(fid,'\n');
# 		fprintf(fid,'%s \n', ['Data file: ' inputFile]); fprintf(fid,'\n');
# 		fprintf(fid,'%s \n', ['Number of clustered individuals: ' ownNum2Str(ninds)]); fprintf(fid,'\n');
# 		fprintf(fid,'%s \n', ['Number of groups in optimal partition: ' ownNum2Str(npops)]); fprintf(fid,'\n');
# 		fprintf(fid,'%s \n', ['Log(marginal likelihood) of optimal partition: ' ownNum2Str(logml)]); fprintf(fid,'\n');
# 	end

# 	cluster_count = length(unique(PARTITION));
# 	disp(['Best Partition: ']);
# 	if (fid ~= -1)
# 		fprintf(fid,'%s \n',['Best Partition: ']); fprintf(fid,'\n');
# 	end
# 	for m=1:cluster_count
# 		indsInM = find(PARTITION==m);
# 		length_of_beginning = 11 + floor(log10(m));
# 		cluster_size = length(indsInM);

# 		if names
# 			text = ['Cluster ' num2str(m) ': {' char(popnames{indsInM(1)})];
# 			for k = 2:cluster_size
# 				text = [text ', ' char(popnames{indsInM(k)})];
# 			end;
# 		else
# 			text = ['Cluster ' num2str(m) ': {' num2str(indsInM(1))];
# 			for k = 2:cluster_size
# 				text = [text ', ' num2str(indsInM(k))];
# 			end;
# 		end
# 		text = [text '}'];
# 		while length(text)>58
# 			%Take one line and display it.
# 			new_line = takeLine(text,58);
# 			text = text(length(new_line)+1:end);
# 			disp(new_line);
# 			if (fid ~= -1)
# 				fprintf(fid,'%s \n',[new_line]);
# 				fprintf(fid,'\n');
# 			end
# 			if length(text)>0
# 				text = [blanks(length_of_beginning) text];
# 			else
# 				text = [];
# 			end;
# 		end;
# 		if ~isempty(text)
# 			disp(text);
# 			if (fid ~= -1)
# 				fprintf(fid,'%s \n',[text]);
# 				fprintf(fid,'\n');
# 			end
# 		end;
# 	end

# 	if npops > 1
# 		disp(' ');
# 		disp(' ');
# 		disp('Changes in log(marginal likelihood) if indvidual i is moved to group j:');
# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', ['Changes in log(marginal likelihood) if indvidual i is moved to group j:']); fprintf(fid, '\n');
# 		end

# 		if names
# 			nameSizes = zeros(ninds,1);
# 			for i = 1:ninds
# 				nimi = char(popnames{i});
# 				nameSizes(i) = length(nimi);
# 			end
# 			maxSize = max(nameSizes);
# 			maxSize = max(maxSize, 5);
# 			erotus = maxSize - 5;
# 			alku = blanks(erotus);
# 			ekarivi = [alku '  ind' blanks(6+erotus)];
# 		else
# 			ekarivi = '  ind      ';
# 		end
# 		for i = 1:cluster_count
# 			ekarivi = [ekarivi ownNum2Str(i) blanks(8-floor(log10(i)))];
# 		end
# 		disp(ekarivi);
# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [ekarivi]); fprintf(fid, '\n');
# 		end

# 		%ninds = size(data,1)/rowsFromInd;
# 		changesInLogml = LOGDIFF';
# 		for ind = 1:ninds
# 			%[muutokset, diffInCounts] = laskeMuutokset(ind, rowsFromInd, data, ...
# 			%    adjprior, priorTerm);
# 			%changesInLogml(:,ind) = muutokset;
# 			muutokset = changesInLogml(:,ind);

# 			if names
# 				nimi = char(popnames{ind});
# 				rivi = [blanks(maxSize - length(nimi)) nimi ':'];
# 			else
# 				rivi = [blanks(4-floor(log10(ind))) ownNum2Str(ind) ':'];
# 			end
# 			for j = 1:npops
# 				rivi = [rivi '  ' logml2String(omaRound(muutokset(j)))];
# 			end
# 			disp(rivi);
# 			if (fid ~= -1)
# 				fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
# 			end
# 		end

# 		disp(' '); disp(' ');
# 		disp('KL-divergence matrix in PHYLIP format:');

# 		dist_mat = zeros(npops, npops);
# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', [' ']); %fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', ['KL-divergence matrix in PHYLIP format:']); %fprintf(fid, '\n');
# 		end

# 		maxnoalle = size(COUNTS,1);
# 		nloci = size(COUNTS,2);
# 		d = zeros(maxnoalle, nloci, npops);
# 		prior = adjprior;
# 		prior(find(prior==1))=0;
# 		nollia = find(all(prior==0));  %Lokukset, joissa oli havaittu vain yht?alleelia.
# 		prior(1,nollia)=1;
# 		for pop1 = 1:npops
# 			d(:,:,pop1) = (squeeze(COUNTS(:,:,pop1))+prior) ./ repmat(sum(squeeze(COUNTS(:,:,pop1))+prior),maxnoalle,1);
# 			%dist1(pop1) = (squeeze(COUNTS(:,:,pop1))+adjprior) ./ repmat((SUMCOUNTS(pop1,:)+adjprior), maxnoalle, 1);
# 		end
# 	%     ekarivi = blanks(7);
# 	%     for pop = 1:npops
# 	%         ekarivi = [ekarivi num2str(pop) blanks(7-floor(log10(pop)))];
# 	%     end
# 		ekarivi = num2str(npops);
# 		disp(ekarivi);
# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [ekarivi]); %fprintf(fid, '\n');
# 		end

# 		for pop1 = 1:npops
# 			% rivi = [blanks(2-floor(log10(pop1))) num2str(pop1) '  '];
# 			for pop2 = 1:pop1-1
# 				dist1 = d(:,:,pop1); dist2 = d(:,:,pop2);
# 				div12 = sum(sum(dist1.*log2((dist1+10^-10) ./ (dist2+10^-10))))/nloci;
# 				div21 = sum(sum(dist2.*log2((dist2+10^-10) ./ (dist1+10^-10))))/nloci;
# 				div = (div12+div21)/2;
# 				% rivi = [rivi kldiv2str(div) '  '];
# 				dist_mat(pop1,pop2) = div;
# 			end
# 			% disp(rivi);
# 			% if (fid ~= -1)
# 			%    fprintf(fid, '%s \n', [rivi]); fprintf(fid, '\n');
# 			% end
# 		end


# 		dist_mat = dist_mat + dist_mat'; % make it symmetric
# 		for pop1 = 1:npops
# 			rivi = ['Cluster_' num2str(pop1) ' '];
# 			for pop2 = 1:npops
# 				rivi = [rivi kldiv2str(dist_mat(pop1,pop2)) ' '];
# 			end
# 			disp(rivi);
# 			if (fid ~= -1)
# 				fprintf(fid, '%s \n', [rivi]); %fprintf(fid, '\n');
# 			end
# 		end
# 	end

# 	disp(' ');
# 	disp(' ');
# 	disp('List of sizes of 10 best visited partitions and corresponding log(ml) values');

# 	if (fid ~= -1)
# 		fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 		fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 		fprintf(fid, '%s \n', ['List of sizes of 10 best visited partitions and corresponding log(ml) values']); fprintf(fid, '\n');
# 	end

# 	partitionSummary = sortrows(partitionSummary,2);
# 	partitionSummary = partitionSummary(size(partitionSummary,1):-1:1 , :);
# 	partitionSummary = partitionSummary(find(partitionSummary(:,2)>-1e49),:);
# 	if size(partitionSummary,1)>10
# 		vikaPartitio = 10;
# 	else
# 		vikaPartitio = size(partitionSummary,1);
# 	end
# 	for part = 1:vikaPartitio
# 		line = [num2str(partitionSummary(part,1)) '    ' num2str(partitionSummary(part,2))];
# 		disp(line);
# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
# 		end
# 	end

# 	if ~fixedK
# 		disp(' ');
# 		disp(' ');
# 		disp('Probabilities for number of clusters');

# 		if (fid ~= -1)
# 			fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', [' ']); fprintf(fid, '\n');
# 			fprintf(fid, '%s \n', ['Probabilities for number of clusters']); fprintf(fid, '\n');
# 		end

# 		npopsTaulu = unique(partitionSummary(:,1));
# 		len = length(npopsTaulu);
# 		probs = zeros(len,1);
# 		partitionSummary(:,2) = partitionSummary(:,2)-max(partitionSummary(:,2));
# 		sumtn = sum(exp(partitionSummary(:,2)));
# 		for i=1:len
# 			npopstn = sum(exp(partitionSummary(find(partitionSummary(:,1)==npopsTaulu(i)),2)));
# 			probs(i) = npopstn / sumtn;
# 		end
# 		for i=1:len
# 			if probs(i)>1e-5
# 				line = [num2str(npopsTaulu(i)) '   ' num2str(probs(i))];
# 				disp(line);
# 				if (fid ~= -1)
# 					fprintf(fid, '%s \n', [line]); fprintf(fid, '\n');
# 				end
# 			end
# 		end
# 	end

# 	if (fid ~= -1)
# 		fclose(fid);
# 	else
# 		diary off
# 	end


# 	%---------------------------------------------------------------


# 	function dispLine;
# 	disp('---------------------------------------------------');

# 	%--------------------------------------------------------------

# 	function num2 = omaRound(num)
# 	% Py�rist�� luvun num 1 desimaalin tarkkuuteen
# 	num = num*10;
# 	num = round(num);
# 	num2 = num/10;

# 	%---------------------------------------------------------


# 	function digit = palautaYks(num,yks)
# 	% palauttaa luvun num 10^yks termin kertoimen
# 	% string:in?
# 	% yks t�ytyy olla kokonaisluku, joka on
# 	% v�hint��n -1:n suuruinen. Pienemmill?
# 	% luvuilla tapahtuu jokin py�ristysvirhe.

# 	if yks>=0
# 		digit = rem(num, 10^(yks+1));
# 		digit = floor(digit/(10^yks));
# 	else
# 		digit = num*10;
# 		digit = floor(rem(digit,10));
# 	end
# 	digit = num2str(digit);


# 	function mjono = kldiv2str(div)
# 	mjono = '      ';
# 	if abs(div)<100
# 		%Ei tarvita e-muotoa
# 		mjono(6) = num2str(rem(floor(div*1000),10));
# 		mjono(5) = num2str(rem(floor(div*100),10));
# 		mjono(4) = num2str(rem(floor(div*10),10));
# 		mjono(3) = '.';
# 		mjono(2) = num2str(rem(floor(div),10));
# 		arvo = rem(floor(div/10),10);
# 		if arvo>0
# 			mjono(1) = num2str(arvo);
# 		end

# 	else
# 		suurinYks = floor(log10(div));
# 		mjono(6) = num2str(suurinYks);
# 		mjono(5) = 'e';
# 		mjono(4) = palautaYks(abs(div),suurinYks-1);
# 		mjono(3) = '.';
# 		mjono(2) = palautaYks(abs(div),suurinYks);
# 	end

# 	%--------------------------------------------------------------------


# 	function newline = takeLine(description,width)
# 	%Returns one line from the description: line ends to the first
# 	%space after width:th mark.
# 	newLine = description(1:width);
# 	n = width+1;
# 	while ~isspace(description(n)) & n<length(description)
# 		n = n+1;
# 	end;
# 	newline = description(1:n);

# 	%-------------------------------------------------------

# 	function pal = testaaPop(rivi)
# 	% pal=1, mik�li rivi alkaa jollain seuraavista
# 	% kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
# 	% tapauksissa pal=0.

# 	if length(rivi)<3
# 		pal = 0;
# 		return
# 	end
# 	if (all(rivi(1:3)=='Pop') | ...
# 		all(rivi(1:3)=='pop') | ...
# 		all(rivi(1:3)=='POP'))
# 		pal = 1;
# 		return
# 	else
# 		pal = 0;
# 		return
# 	end



# 	%----------------------------------------------------------------------------


# 	function dist2 = laskeOsaDist(inds2, dist, ninds)
# 	% Muodostaa dist vektorista osavektorin, joka sis�lt�� yksil�iden inds2
# 	% v�liset et�isyydet. ninds=kaikkien yksil�iden lukum��r?

# 	ninds2 = length(inds2);
# 	apu = zeros(nchoosek(ninds2,2),2);
# 	rivi = 1;
# 	for i=1:ninds2-1
# 		for j=i+1:ninds2
# 			apu(rivi, 1) = inds2(i);
# 			apu(rivi, 2) = inds2(j);
# 			rivi = rivi+1;
# 		end
# 	end
# 	apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
# 	dist2 = dist(apu);

# 	%--------------------------------------------------------------------------

# 	function [emptyPop, pops] = findEmptyPop(npops)
# 	% Palauttaa ensimm�isen tyhj�n populaation indeksin. Jos tyhji?
# 	% populaatioita ei ole, palauttaa -1:n.

# 	global PARTITION;
# 	pops = unique(PARTITION)';
# 	if (length(pops) ==npops)
# 		emptyPop = -1;
# 	else
# 		popDiff = diff([0 pops npops+1]);
# 		emptyPop = min(find(popDiff > 1));
# 	end

# 	%------------------------------------------------------


# 	function [data, popnames] = lueGenePopData(tiedostonNimi)

# 	fid = fopen(tiedostonNimi);
# 	line = fgetl(fid);  %ensimm�inen rivi
# 	line = fgetl(fid);  %toinen rivi
# 	count = rivinSisaltamienMjonojenLkm(line);

# 	line = fgetl(fid);
# 	lokusRiveja = 1;
# 	while (testaaPop(line)==0)
# 		lokusRiveja = lokusRiveja+1;
# 		line = fgetl(fid);
# 	end

# 	if lokusRiveja>1
# 		nloci = lokusRiveja;
# 	else
# 		nloci = count;
# 	end

# 	popnames = cell(10,2);
# 	data = zeros(100, nloci+1);
# 	nimienLkm=0;
# 	ninds=0;
# 	poimiNimi=1;
# 	digitFormat = -1;
# 	while line ~= -1
# 		line = fgetl(fid);

# 		if poimiNimi==1
# 			%Edellinen rivi oli 'pop'
# 			nimienLkm = nimienLkm+1;
# 			ninds = ninds+1;
# 			if nimienLkm>size(popnames,1);
# 				popnames = [popnames; cell(10,2)];
# 			end
# 			nimi = lueNimi(line);
# 			if digitFormat == -1
# 				digitFormat = selvitaDigitFormat(line);
# 				divider = 10^digitFormat;
# 			end
# 			popnames{nimienLkm, 1} = {nimi};   %N�in se on greedyMix:iss�kin?!?
# 			popnames{nimienLkm, 2} = ninds;
# 			poimiNimi=0;

# 			data = addAlleles(data, ninds, line, divider);

# 		elseif testaaPop(line)
# 			poimiNimi = 1;

# 		elseif line ~= -1
# 			ninds = ninds+1;
# 			data = addAlleles(data, ninds, line, divider);
# 		end
# 	end

# 	data = data(1:ninds*2,:);
# 	popnames = popnames(1:nimienLkm,:);
# 	fclose(fid);

# 	%--------------------------------------------------------


# 	function data = addAlleles(data, ind, line, divider)
# 	% Lisaa BAPS-formaatissa olevaan datataulukkoon
# 	% yksil�� ind vastaavat rivit. Yksil�n alleelit
# 	% luetaan genepop-formaatissa olevasta rivist?
# 	% line. Jos data on 3 digit formaatissa on divider=1000.
# 	% Jos data on 2 digit formaatissa on divider=100.

# 	nloci = size(data,2)-1;
# 	if size(data,1) < 2*ind
# 		data = [data; zeros(100,nloci+1)];
# 	end

# 	k=1;
# 	merkki=line(k);
# 	while ~isequal(merkki,',')
# 	   k=k+1;
# 	   merkki=line(k);
# 	end
# 	line = line(k+1:end);
# 	clear k; clear merkki;

# 	alleeliTaulu = sscanf(line,'%d');

# 	if length(alleeliTaulu)~=nloci
# 		disp('Incorrect data format.');
# 	end

# 	for j=1:nloci
# 		ekaAlleeli = floor(alleeliTaulu(j)/divider);
# 		if ekaAlleeli==0 ekaAlleeli=-999; end;
# 		tokaAlleeli = rem(alleeliTaulu(j),divider);
# 		if tokaAlleeli==0 tokaAlleeli=-999; end

# 		data(2*ind-1,j) = ekaAlleeli;
# 		data(2*ind,j) = tokaAlleeli;
# 	end

# 	data(2*ind-1,end) = ind;
# 	data(2*ind,end) = ind;

# 	%------------------------------------------------------


# 	function count = rivinSisaltamienMjonojenLkm(line)
# 	% Palauttaa line:n sis�lt�mien mjonojen lukum��r�n.
# 	% Mjonojen v�liss?t�ytyy olla v�lily�nti.
# 	count = 0;
# 	pit = length(line);
# 	tila = 0;    %0, jos odotetaan v�lily�ntej? 1 jos odotetaan muita merkkej?
# 	for i=1:pit
# 		merkki = line(i);
# 		if (isspace(merkki) & tila==0)
# 			%Ei tehd?mit��n.
# 		elseif (isspace(merkki) & tila==1)
# 			tila = 0;
# 		elseif (~isspace(merkki) & tila==0)
# 			tila = 1;
# 			count = count+1;
# 		elseif (~isspace(merkki) & tila==1)
# 			%Ei tehd?mit��n
# 		end
# 	end


# 	%-------------------------------------------------------

# 	function nimi = lueNimi(line)
# 	%Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
# 	n = 1;
# 	merkki = line(n);
# 	nimi = '';
# 	while ~isequal(merkki,',')
# 		nimi = [nimi merkki];
# 		n = n+1;
# 		merkki = line(n);
# 	end

# 	%-------------------------------------------------------

# 	function df = selvitaDigitFormat(line)
# 	% line on ensimm�inen pop-sanan j�lkeinen rivi
# 	% Genepop-formaatissa olevasta datasta. funktio selvitt��
# 	% rivin muodon perusteella, ovatko datan alleelit annettu
# 	% 2 vai 3 numeron avulla.

# 	n = 1;
# 	merkki = line(n);
# 	while ~isequal(merkki,',')
# 		n = n+1;
# 		merkki = line(n);
# 	end

# 	while ~any(merkki == '0123456789');
# 		n = n+1;
# 		merkki = line(n);
# 	end
# 	numeroja = 0;
# 	while any(merkki == '0123456789');
# 		numeroja = numeroja+1;
# 		n = n+1;
# 		merkki = line(n);
# 	end

# 	df = numeroja/2;



# 	function loggis = laskeLoggis(counts, sumcounts, adjprior)
# 	npops = size(counts,3);

# 	logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
# 		- npops*sum(sum(gammaln(adjprior))) - ...
# 		sum(sum(gammaln(1+sumcounts)));
# 	loggis = logml2;