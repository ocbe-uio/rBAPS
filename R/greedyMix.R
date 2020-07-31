#' @title Clustering of individuals
#' @param tietue File
#' @param format Format of the data ("BAPS", "GenePop" or "Preprocessed")
#' @param savePreProcessed Save the pre-processed data?
#' @param filePreProcessed Is the file already processed?
#' @importFrom utils read.delim
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
		if (length(input_type_name) == 0) {
			stop('Invalid alternative')
		} else if (input_type_name == 'BAPS-format') {
			# ------------------------------------------------------------------
			# Treating BAPS-formatted files
			# ------------------------------------------------------------------
			if (!is(tietue, "character")) {
				pathname_filename <- uigetfile(
					"*.txt", "Load data in BAPS-format"
				)
			} else {
				pathname_filename <- tietue
			}

			# ASK: remove?
			# if ~isempty(partitionCompare)
			# 	fprintf(1,'Data: %s\n',[pathname filename]);
			# end

			data <- read.delim(pathname_filename) # TODO: discover delimiter
			ninds <- testaaOnkoKunnollinenBapsData(data)  # testing
			if (ninds == 0) stop('Incorrect Data-file')

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
		} else if (input_type_name == 'GenePop-format') {
			# ------------------------------------------------------------------
			# Treating GenePop-formatted files
			# ------------------------------------------------------------------
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

			kunnossa <- testaaGenePopData(filename_pathname)
			if (kunnossa == 0) stop("testaaGenePopData returned 0")
			data_popnames <- lueGenePopData(filename_pathname)
			data          <- data_popnames$data
			popnames      <- data_popnames$popnames

# 			h0 = findobj('Tag','filename1_text');
# 			set(h0,'String',filename); clear h0;

			list_dranap <- handleData(data)
			data        <- list_dranap$newData
			rowsFromInd <- list_dranap$rowsFromInd
			alleleCodes <- list_dranap$alleleCodes
			noalle      <- list_dranap$noalle
			adjprior    <- list_dranap$adjprior
			priorTerm   <- list_dranap$prioterm

			list_Zd <- newGetDistances(data,rowsFromInd) # FIXME: debug
			Z       <- list_Zd$Z
			dist    <- list_Zd$dist

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
		} else if (input_type_name == 'Preprocessed data') {
			# ------------------------------------------------------------------
			# Handling Pre-processed data
			# ------------------------------------------------------------------
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
	c             <- list()
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
	# if (!is.null(partitionCompare)) {
	# 	nsamplingunits <- size(c$rows, 1)
	# 	partitions     <- partitionCompare$partitions
	# 	npartitions    <- size(partitions, 2)
	# 	partitionLogml <- zeros(1, npartitions)
	# 	for (i in seq_len(npartitions)) {
	# 		# number of unique partition lables
	# 		npops <- length(unique(partitions[, i]))

	# 		partitionInd <- zeros(ninds * rowsFromInd, 1)
	# 		partitionSample <- partitions[, i]
	# 		for (j in seq_len(nsamplingunits)) {
	# 			partitionInd[c$rows[j, 1]:c$rows[j, 2]] <- partitionSample[j]
	# 		}
	# 		# partitionLogml[i] = initialCounts(
	# 		# 	partitionInd,
	# 		# 	data[, seq_len(end - 1)],
	# 		# 	npops,
	# 		# 	c$rows,
	# 		# 	noalle,
	# 		# 	adjprior
	# 		# ) #TODO translate
	# 	}
	# 	# return the logml result
	# 	partitionCompare$logmls <- partitionLogml
	# 	# set(h1, 'userdata', partitionCompare) # ASK remove?
	# 	return()
	# }

	# ASK remove (graphical part)?
	# if (fixedK) {
	# 	#logml_npops_partitionSummary <- indMix_fixK(c) # ASK translate?
	# } else {
	# 	#logml_npops_partitionSummary <- indMix(c) # ASK translate?
	# }
	# if (logml_npops_partitionSummary$logml == 1) return()

	data <- data[, seq_len(ncol(data) - 1)]

	# ASK: remove?
	# h0 = findobj('Tag','filename1_text')
	# inp = get(h0,'String');
	# h0 = findobj('Tag','filename2_text')
	# outp = get(h0,'String');

	browser() # TEMP
	changesInLogml <- writeMixtureInfo(
		logml, rowsFromInd, data, adjprior, priorTerm, outp, inp,
		popnames, fixedK
	) # FIXME: broken

	# viewMixPartition(PARTITION, popnames) # ASK translate? On graph folder

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

		# ======================================================================
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
		# ======================================================================

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

# 	%---------------------------------------------------------------

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