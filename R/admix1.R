#' @title Admixture analysis
#' @param tietue a named record list
#' @details If the record == -1, the mixture results file is loaded. Otherwise,
#' will the required variables be retrieved from the record fields?
#' `tietue`should contain the following elements: PARTITION, COUNTS, SUMCOUNTS,
#' alleleCodes, adjprior, popnames, rowsFromInd, data, npops, noalle
#' @param tietue tietue
#' @importFrom methods is
#' @export
admix1 <- function(tietue) {
    if (!is.list(tietue)) {
        message('Load mixture result file. These are the files in this directory:')
        print(list.files())
        pathname_filename <- file.choose()
        if (!file.exists(pathname_filename)) {
            stop(
                "File ", pathname_filename,
                " does not exist. Check spelling and location."
            )
        } else {
            cat('---------------------------------------------------\n');
            message('Reading mixture result from: ', pathname_filename, '...')
        }
        Sys.sleep(0.0001) #TODO: remove

        # ASK: what is this supposed to do? What do graphic obj have to do here?
        # h0 = findobj('Tag','filename1_text');
        # set(h0,'String',filename); clear h0;

        struct_array <- load(pathname_filename)
        if (isfield(struct_array, 'c')) { #Matlab versio
            c <- struct_array$c
            if (!isfield(c, 'PARTITION') | !isfield(c,'rowsFromInd')) {
                stop('Incorrect file format')
            }
        } else if (isfield(struct_array, 'PARTITION')) {  #Mideva versio
            c <- struct_array
            if (!isfield(c,'rowsFromInd')) stop('Incorrect file format')
        } else {
            stop('Incorrect file format')
        }

        if (isfield(c, 'gene_lengths') &
            strcmp(c$mixtureType, 'linear_mix') |
            strcmp(c$mixtureType, 'codon_mix')) { # if the mixture is from a linkage model
            # Redirect the call to the linkage admixture function.
            # call function noindex to remove the index column
            c$data <- noIndex(c$data, c$noalle)
            # linkage_admix(c) # TODO: obsolete. remove.
            # return
            stop("linkage_admix not implemented")
        }
        PARTITION <- c$PARTITION
        COUNTS <- c$COUNTS
        SUMCOUNTS <- c$SUMCOUNTS
        alleleCodes <- c$alleleCodes
        adjprior <- c$adjprior
        popnames <- c$popnames
        rowsFromInd <- c$rowsFromInd
        data <- c$data
        npops <- c$npops
        noalle <- c$noalle
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

    answers <- inputdlg(
        prompt = paste(
            "Input the minimum size of a population that will",
            "be taken into account when admixture is estimated."
        ),
        definput = 5
    )
    alaRaja <- as.numeric(answers)
    npops <- poistaLiianPienet(npops, rowsFromInd, alaRaja)

    nloci <- size(COUNTS, 2)
    ninds <- size(data, 1) / rowsFromInd

    answers <- inputdlg('Input number of iterations', definput=50)
    if (isempty(answers)) return()
    iterationCount <- as.numeric(answers[1, 1]) # maybe [[]]?

    answers <- inputdlg(
        prompt = 'Input number of reference individuals from each population',
        definput = 50
    )
    if (isempty(answers)) {
        nrefIndsInPop <- 50
    } else {
        nrefIndsInPop <- as.numeric(answers[1, 1])
    }

    answers <- inputdlg(
        prompt = 'Input number of iterations for reference individuals',
        definput = 10
    )
    if (isempty(answers)) return()
    iterationCountRef <- as.numeric(answers[1, 1])

    # First calculate log-likelihood ratio for all individuals:
    likelihood <- zeros(ninds, 1)
    allfreqs <- computeAllFreqs2(noalle)
    for (ind in 1:ninds) {
        omaFreqs <- computePersonalAllFreqs(ind, data, allfreqs, rowsFromInd)
        osuusTaulu <- zeros(1, npops)
        if (PARTITION[ind] == 0) {
            # Yksil?on outlier
        } else if (PARTITION[ind] != 0) {
            if (PARTITION[ind] > 0) {
                osuusTaulu[PARTITION[ind]] <- 1
            } else {
                # Yksilöt, joita ei ole sijoitettu mihinkään koriin.
                arvot <- zeros(1, npops)
                for (q in 1:npops) {
                    osuusTaulu <- zeros(1, npops)
                    osuusTaulu[q] <- 1
                    arvot[q] <- computeIndLogml(omaFreqs, osuusTaulu)
                }
                iso_arvo <- max(arvot)
                isoimman_indeksi <- match(max(arvot), arvot)
                osuusTaulu <- zeros(1, npops)
                osuusTaulu[isoimman_indeksi] <- 1
                PARTITION[ind] <- isoimman_indeksi
            }
            logml <- computeIndLogml(omaFreqs, osuusTaulu)
            logmlAlku <- logml
            for (osuus in c(0.5, 0.25, 0.05, 0.01)) {
                etsiResult <- etsiParas(osuus, osuusTaulu, omaFreqs, logml)
                osuusTaulu <- etsiResult[1]
                logml <- etsiResult[2]
            }
            logmlLoppu <- logml
            likelihood[ind] <- logmlLoppu - logmlAlku
        }
    }

    # Analyze further only individuals who have log-likelihood ratio larger than 3:
    to_investigate <- t(find(likelihood > 3))
    cat('Possibly admixed individuals:\n')
    for (i in 1:length(to_investigate)) {
        cat(as.character(to_investigate[i]))
    }
    cat(' ')
    cat('Populations for possibly admixed individuals:\n')
    admix_populaatiot <- unique(PARTITION[to_investigate])
    for (i in 1:length(admix_populaatiot)) {
        cat(as.character(admix_populaatiot[i]))
    }

    # THUS, there are two types of individuals, who will not be analyzed with
    # simulated allele frequencies: those who belonged to a mini-population
    # which was removed, and those who have log-likelihood ratio less than 3.
    # The value in the PARTITION for the first kind of individuals is 0. The
    # second kind of individuals can be identified, because they do not
    # belong to "to_investigate" array. When the results are presented, the
    # first kind of individuals are omitted completely, while the second kind
    # of individuals are completely put to the population, where they ended up
    # in the mixture analysis. These second type of individuals will have a
    # unit p-value.


    # Simulate allele frequencies a given number of times and save the average
    # result to "proportionsIt" array.

    proportionsIt <- zeros(ninds, npops)
    for (iterationNum in 1:iterationCount) {
        cat('Iter:', as.character(iterationNum))
        allfreqs <- simulateAllFreqs(noalle)  # Allele frequencies on this iteration.

        for (ind in to_investigate) {
            #disp(num2str(ind));
            omaFreqs <- computePersonalAllFreqs(
                ind, data, allfreqs, rowsFromInd
            )
            osuusTaulu <- zeros(1, npops)
            if (PARTITION[ind] == 0) {
                # Yksil?on outlier
            } else if (PARTITION[ind] != 0) {
                if (PARTITION[ind] > 0) {
                    osuusTaulu[PARTITION[ind]] <- 1
                } else {
                    # Yksilöt, joita ei ole sijoitettu mihinkään koriin.
                    arvot <- zeros(1, npops)
                    for (q in 1:npops) {
                        osuusTaulu <- zeros(1, npops)
                        osuusTaulu[q] <- 1
                        arvot[q] <- computeIndLogml(omaFreqs, osuusTaulu)
                    }
                    iso_arvo <- max(arvot)
                    isoimman_indeksi <- match(max(arvot), arvot)
                    osuusTaulu <- zeros(1, npops)
                    osuusTaulu[isoimman_indeksi] <- 1
                    PARTITION[ind] <- isoimman_indeksi
                }
                logml <- computeIndLogml(omaFreqs, osuusTaulu)

                for (osuus in c(0.5, 0.25, 0.05, 0.01)) {
                    etsiResult <- etsiParas(osuus, osuusTaulu, omaFreqs, logml)
                    osuusTaulu <- etsiResult[1]
                    logml <- etsiResult[2]
                }
            }
            proportionsIt[ind, ] <- proportionsIt[ind, ] * (iterationNum - 1) +
                osuusTaulu
            proportionsIt[ind, ] <- proportionsIt[ind, ] / iterationNum
        }
    }

    #disp(['Creating ' num2str(nrefIndsInPop) ' reference individuals from ']);
    #disp('each population.');

    #allfreqs = simulateAllFreqs(noalle);  # Simuloidaan alleelifrekvenssisetti
    allfreqs <- computeAllFreqs2(noalle); # Koitetaan tällaista.


    # Initialize the data structures, which are required in taking the missing
    # data into account:
    n_missing_levels <- zeros(npops, 1)  # number of different levels of "missingness" in each pop (max 3).
    missing_levels <- zeros(npops, 3)  # the mean values for different levels.
    missing_level_partition <- zeros(ninds, 1) # level of each individual (one of the levels of its population).
    for (i in 1:npops) {
        inds <- find(PARTITION == i)
        # Proportions of non-missing data for the individuals:
        non_missing_data <- zeros(length(inds), 1)
        for (j in 1:length(inds)) {
            ind <- inds[j]
            non_missing_data[j] <- length(
                find(data[(ind - 1) * rowsFromInd + 1:ind * rowsFromInd, ] > 0)
            ) / (rowsFromInd * nloci)
        }
        if (all(non_missing_data > 0.9)) {
            n_missing_levels[i] <- 1
            missing_levels[i, 1] <- mean(non_missing_data)
            missing_level_partition[inds] <- 1
        } else {
            # TODO: fix syntax
            # [ordered, ordering] = sort(non_missing_data);
            ordered <- ordering <- sort(non_missing_data)
            #part = learn_simple_partition(ordered, 0.05);
            part <- learn_partition_modified(ordered)
            aux <- sortrows(cbind(part, ordering), 2)
            part = aux[, 1]
            missing_level_partition[inds]<- part
            n_levels <- length(unique(part))
            n_missing_levels[i] <- n_levels
            for (j in 1:n_levels) {
                missing_levels[i, j] <- mean(non_missing_data[find(part == j)])
            }
        }
    }

    # Create and analyse reference individuals for populations
    # with potentially admixed individuals:
    refTaulu <- zeros(npops, 100, 3)
    for (pop in t(admix_populaatiot)) {

        for (level in 1:n_missing_levels[pop]) {

            potential_inds_in_this_pop_and_level <-
                find(
                    PARTITION == pop & missing_level_partition == level &
                    likelihood > 3
                ) # Potential admix individuals here.

            if (!isempty(potential_inds_in_this_pop_and_level)) {

                #refData = simulateIndividuals(nrefIndsInPop,rowsFromInd,allfreqs);
                refData <- simulateIndividuals(
                    nrefIndsInPop, rowsFromInd, allfreqs, pop,
                    missing_levels[pop, level]
                )

                cat(
                    'Analysing the reference individuals from pop', pop,
                    '(level', level, ').'
                )
                refProportions <- zeros(nrefIndsInPop, npops)
                for (iter in 1:iterationCountRef) {
                    #disp(['Iter: ' num2str(iter)]);
                    allfreqs <- simulateAllFreqs(noalle)

                    for (ind in 1:nrefIndsInPop) {
                        omaFreqs <- computePersonalAllFreqs(
                            ind, refData, allfreqs, rowsFromInd
                        )
                        osuusTaulu <- zeros(1, npops)
                        osuusTaulu[pop] <- 1
                        logml <- computeIndLogml(omaFreqs, osuusTaulu)
                        for (osuus in c(0.5, 0.25, 0.05, 0.01)) {
                            etsiResult <- etsiParas(
                                osuus, osuusTaulu, omaFreqs, logml
                            )
                            osuusTaulu <- etsiResult[1]
                            logml <- etsiResult[2]
                        }
                        refProportions[ind, ] <-
                            refProportions[ind, ] * (iter - 1) + osuusTaulu
                        refProportions[ind, ] <- refProportions[ind, ] / iter
                        }
                }
                for (ind in 1:nrefIndsInPop) {
                    omanOsuus <- refProportions[ind, pop]
                    if (round(omanOsuus * 100) == 0) {
                        omanOsuus <- 0.01
                    }
                    if (abs(omanOsuus) < 1e-5) {
                        omanOsuus <- 0.01
                    }
                    refTaulu[pop, round(omanOsuus*100), level] <-
                        refTaulu[pop, round(omanOsuus*100),level] + 1
                }
            }
        }
    }

    # Rounding of the results:
    proportionsIt <- proportionsIt * 100
    proportionsIt <- round(proportionsIt)
    proportionsIt <- proportionsIt / 100
    for (ind in 1:ninds) {
        if (!any(to_investigate == ind)) {
            if (PARTITION[ind] > 0) {
                proportionsIt[ind, PARTITION[ind]] <- 1
            }
        } else {
            # In case of a rounding error, the sum is made equal to unity by
            # fixing the largest value.
            if ((PARTITION[ind] > 0) & (sum(proportionsIt[ind, ]) != 1)) {
                isoin <- max(proportionsIt[ind, ])
                indeksi <- match(isoin, max(proportionsIt[ind, ]))
                erotus <- sum(proportionsIt[ind, ]) - 1
                proportionsIt[ind, indeksi] <- isoin - erotus
            }
        }
    }

    # Calculate p-value for each individual:
    uskottavuus <- zeros(ninds, 1)
    for (ind in 1:ninds) {
        pop <- PARTITION[ind]
        if (pop == 0) { # Individual is outlier
            uskottavuus[ind] <- 1
        } else if (isempty(find(to_investigate == ind))) {
            # Individual had log-likelihood ratio<3
            uskottavuus[ind] <- 1
        } else {
            omanOsuus <- proportionsIt[ind, pop]
            if (abs(omanOsuus) < 1e-5) {
                omanOsuus <- 0.01
            }
            if (round(omanOsuus*100)==0) {
                omanOsuus <- 0.01
            }
            level <- missing_level_partition[ind]
            refPienempia <- sum(refTaulu[pop, 1:round(100*omanOsuus), level])
            uskottavuus[ind] <- refPienempia / nrefIndsInPop
        }
    }

    tulostaAdmixtureTiedot(proportionsIt, uskottavuus, alaRaja, iterationCount) # TODO: textual outputs. probably not necessary. translate nonetheless
    viewPartition(proportionsIt, popnames)  # TODO: adapt

    talle = inputdlg('Do you want to save the admixture results? [Y/n]', 'y')
    if (talle %in% c('y', 'Y', 'yes', 'Yes')) {
        #waitALittle;
        filename <- inputdlg(
            'Save results as (file name):', 'admixture_results.rda'
        )


        if (filename == 0) {
            # Cancel was pressed
            return()
        } else { # copy 'baps4_output.baps' into the text file with the same name.
            if (file.exists('baps4_output.baps')) {
                file.copy('baps4_output.baps', paste0(filename, '.txt'))
                file.remove('baps4_output.baps')
            }
        }

        if (!is(tietue, "list")) {
            c$proportionsIt <- proportionsIt
            c$pvalue <- uskottavuus # Added by Jing
            c$mixtureType <- 'admix' # Jing
            c$admixnpops <- npops;
            save(c, file=filename)
        } else {
            tietue$proportionsIt <- proportionsIt
            tietue$pvalue <- uskottavuus; # Added by Jing
            tietue$mixtureType <- 'admix'
            tietue$admixnpops <- npops
            save(tietue, file=filename)
        }
    }
}