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
	# ASK: Unclear when fixedk == TRUE. Remove?
	# check whether fixed k mode is selected
	# h0 <- findobj('Tag','fixk_menu')
	# fixedK = get(h0, 'userdata');
	fixedK <- FALSE

	# if fixedK
	# 	if ~(fixKWarning == 1) % call function fixKWarning
	# 		return
	# 	end
	# end

	# ASK: ditto
	# % check whether partition compare mode is selected
	# h1 = findobj('Tag','partitioncompare_menu');
	# partitionCompare = get(h1, 'userdata');
	partitionCompare <- FALSE

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

			data <- read.delim(pathname_filename, header = FALSE, sep = " ")
			data <- as.matrix(data)
			ninds <- testaaOnkoKunnollinenBapsData(data)  # testing
			if (ninds == 0) stop('Incorrect Data-file')

			# ASK: remove?
			# h0 = findobj('Tag','filename1_text');
			# set(h0,'String',filename); clear h0;
			message(
				'When using data which are in BAPS-format, ',
				'you can specify the sampling populations of the ',
				'individuals by giving two additional files: ',
				'one containing the names of the populations, ',
				'the other containing the indices of the first ',
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

			temp_handleData <- handleData(data)
			data <- temp_handleData$newData
			rowsFromInd <- temp_handleData$rowsFromInd
			alleleCodes <- temp_handleData$alleleCodes
			noalle <- temp_handleData$noalle
			adjprior <- temp_handleData$adjprior
			priorTerm <- temp_handleData$priorTerm
			Z_dist <- newGetDistances(data,rowsFromInd)
			Z <- Z_dist$Z
			dist <- Z_dist$dist
			rm(temp_handleData, Z_dist)
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
	# Declaring global variables and changing environment of children functions
	# ==========================================================================
	PARTITION       <- vector()
	COUNTS          <- vector()
	SUMCOUNTS       <- vector()
	POP_LOGML       <- vector()
	clearGlobalVars()

	environment(writeMixtureInfo) <- environment()
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
	c$rows <- t(rbind(ekat, ekat + rowsFromInd - 1))

	# ASK remove?
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
	# 	# set(h1, 'userdata', partitionCompare)
	# 	return()
	# }

	if (fixedK) {
		# logml_npops_partitionSummary <- indMix_fixK(c) # TODO: translate
		# logml <- logml_npops_partitionSummary$logml
		# npops <- logml_npops_partitionSummary$npops
		# partitionSummary <- logml_npops_partitionSummary$partitionSummary
	} else {
		logml_npops_partitionSummary <- indMix(c)
		logml <- logml_npops_partitionSummary$logml
		npops <- logml_npops_partitionSummary$npops
		partitionSummary <- logml_npops_partitionSummary$partitionSummary
	}
	if (logml_npops_partitionSummary$logml == 1) return()

	data <- data[, seq_len(ncol(data) - 1)]

	# ASK: remove?
	# h0 = findobj('Tag','filename1_text')
	# inp = get(h0,'String');
	# h0 = findobj('Tag','filename2_text')
	# outp = get(h0,'String');
	inp <- vector()
	outp <- vector()

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
