indMix <- function(c, npops, dispText) {
	# Greedy search algorithm with unknown number of classes for regular
	# clustering.
	# Input npops is not used if called by greedyMix or greedyPopMix.

	logml <- 1
	clearGlobalVars()

	noalle <- c$noalle
	rows <- c$rows
	data <- c$data

	adjprior <- c$adjprior
	priorTerm <- c$priorTerm
	rowsFromInd <- c$rowsFromInd

	if (isfield(c, 'dist')) {
		dist <- c$dist
		Z <- c$Z
	}

	rm(c)
	nargin <- length(as.list(match.call())) - 1

	if (nargin < 2) {
		dispText <- 1
		npopstext <- matrix()
		ready <- FALSE
		teksti <- 'Input upper bound to the number of populations (possibly multiple values)'
		while (!ready) {
			npopstextExtra <- inputdlg(
				teksti,
				1,
				'20'
			)
			if (isempty(npopstextExtra)) { # Painettu Cancel:ia
				return()
			}
			npopstextExtra <- npopstextExtra[1]
			if (length(npopstextExtra)>=255) {
				npopstextExtra <- npopstextExtra[1:255]
				npopstext <- c(npopstext, ' ', npopstextExtra)
				teksti <- 'The input field length limit (255 characters) was reached. Input more values: '
			} else {
				npopstext <- c(npopstext, ' ', npopstextExtra)
				ready <- TRUE
			}
		}
		rm(ready, teksti)
		if (isempty(npopstext) | length(npopstext) == 1) {
			return()
		} else {
			npopsTaulu <- as.numeric(npopstext)
			ykkoset <- find(npopsTaulu == 1)
			npopsTaulu(ykkoset) <- list() # Mik�li ykk�si� annettu yl�rajaksi, ne poistetaan.
			if (isempty(npopsTaulu)) {
				logml <- 1
				partitionSummary <- 1
				npops <- 1
				return()
			}
			rm(ykkoset)
		}
	} else {
		npopsTaulu <- npops
	}

	nruns <- length(npopsTaulu)

	initData <- data
	data <- data[,1:(end - 1)]

	logmlBest <- -1e50
	partitionSummary <- -1e50 * ones(30, 2) # Tiedot 30 parhaasta partitiosta (npops ja logml)
	partitionSummary[,1] <- zeros(30, 1)
	worstLogml <- -1e50
	worstIndex <- 1
	for (run in 1:nruns) {
		npops <- npopsTaulu(run)
		if (dispText) {
			dispLine()
			print(
				paste0(
					'Run ', num2str(run), '/', num2str(nruns),
					', maximum number of populations ', num2str(npops), '.'
				)
			)
		}
		ninds <- size(rows, 1)

		initialPartition <- admixture_initialization(initData, npops, Z) # TODO: translate
		sumcounts_counts_logml = initialCounts(
			initialPartition, data, npops, rows, noalle, adjprior
		) # TODO: translate
		sumcounts <- sumcounts_counts_logml$sumcounts
		counts <- sumcounts_counts_logml$counts
		logml <- sumcounts_counts_logml$logml

		PARTITION <- zeros(ninds, 1)
		for (i in 1:ninds) {
			apu <- rows[i]
			PARTITION[i] <- initialPartition(apu[1])
		}

		COUNTS <- counts
		SUMCOUNTS <- sumcounts
		POP_LOGML <- computePopulationLogml(1:npops, adjprior, priorTerm) # TODO: translate
		LOGDIFF <- repmat(-Inf, c(ninds, npops))
		rm(initialPartition, counts, sumcounts)

		# PARHAAN MIXTURE-PARTITION ETSIMINEN
		nRoundTypes <- 7
		kokeiltu <- zeros(nRoundTypes, 1)
		roundTypes <- c(1, 1)  # Ykk�svaiheen sykli kahteen kertaan.
		ready <- 0
		vaihe <- 1

		if (dispText) {
			print(' ')
			print(
				paste0(
					'Mixture analysis started with initial',
					num2str(npops),
					'populations.'
				)
			)
		}

		while (ready != 1) {
			muutoksia <- 0

			if (dispText) {
				print(paste('Performing steps:', num2str(roundTypes)))
			}

			for (n in 1:length(roundTypes)) {

				round <- roundTypes[n]
				kivaluku <- 0

				if (kokeiltu(round) == 1) { #Askelta kokeiltu viime muutoksen j�lkeen

				} else if (round == 0 | round == 1) { #Yksil�n siirt�minen toiseen populaatioon.
					inds <- 1:ninds
					aputaulu <- c(t(inds), rand(ninds, 1))
					aputaulu <- sortrows(aputaulu, 2)
					inds <- t(aputaulu[, 1])
					muutosNyt <- 0

					for (ind in inds) {
						i1 <- PARTITION[ind]
						muutokset_diffInCounts = laskeMuutokset(
							ind, rows, data, adjprior, priorTerm
						)
						muutokset <- muutokset_diffInCounts$muutokset
						diffInCounts <- muutokset_diffInCounts$diffInCounts

						if (round == 1) {
							maxMuutos <- max_MATLAB(muutokset)[[1]]
							i2 <- max_MATLAB(muutokset)[[2]]
						}

						if (i1 != i2 & maxMuutos > 1e-5) {
							# Tapahtui muutos
							muutoksia <- 1
							if (muutosNyt == 0) {
								muutosNyt <- 1
								if (dispText) {
									print('Action 1')
								}
							}
							kokeiltu <- zeros(nRoundTypes, 1)
							kivaluku <- kivaluku + 1
							updateGlobalVariables(
								ind, i2, diffInCounts, adjprior, priorTerm
							)
							logml <- logml+maxMuutos
							if (logml > worstLogml) {
								partitionSummary_added = addToSummary(
									logml, partitionSummary, worstIndex
								)
								partitionSummary_added <- partitionSummary_added$partitionSummary
								added <- partitionSummary_added$added
								if (added == 1) {
									worstLogml <- min_MATLAB(partitionSummary[, 2])[[1]]
									worstIndex <- min_MATLAB(partitionSummary[, 2])[[2]]
								}
							}
						}
					}

					if (muutosNyt == 0) {
						kokeiltu[round] <- 1
					}

				} else if (round == 2) { # Populaation yhdist�minen toiseen.
					maxMuutos <- 0
					for (pop in 1:npops) {
						muutokset_diffInCounts <- laskeMuutokset2(
							pop, rows, data, adjprior, priorTerm
						)
						muutokset <- muutokset_diffInCounts$muutokset
						diffInCounts <- muutokset_diffInCounts$diffInCounts
						isoin <- max_MATLAB(muutokset)[[1]]
						indeksi <- max_MATLAB(muutokset)[[2]]
						if (isoin > maxMuutos) {
							maxMuutos <- isoin
							i1 <- pop
							i2 <- indeksi
							diffInCountsBest <- diffInCounts
						}
					}

					if (maxMuutos > 1e-5) {
						muutoksia <- 1
						kokeiltu <- zeros(nRoundTypes, 1)
						updateGlobalVariables2(
							i1, i2, diffInCountsBest, adjprior, priorTerm
						)
						logml <- logml + maxMuutos
						if (dispText) {
							print('Action 2')
						}
						if (logml > worstLogml) {
							partitionSummary_added <- addToSummary(
								logml, partitionSummary, worstIndex
							)
							partitionSummary <- partitionSummary_added$partitionSummary
							added <- partitionSummary_added$added
							if (added==1) {
								worstLogml <- min_MATLAB(partitionSummary[, 2])[[1]]
								worstIndex <- min_MATLAB(partitionSummary[, 2])[[2]]
							}
						}
					} else {
						kokeiltu[round] <- 1
					}


				} else if (round == 3 || round == 4) { #Populaation jakaminen osiin.
					maxMuutos <- 0
					ninds <- size(rows, 1)
					for (pop in 1:npops) {
						inds2 <- find(PARTITION == pop)
						ninds2 <- length(inds2)
						if (ninds2 > 2) {
							dist2 <- laskeOsaDist(inds2, dist, ninds)
							Z2 <- linkage(t(dist2))
							if (round == 3) {
								npops2 <- max(min(20, floor(ninds2 / 5)), 2)
							} else if (round == 4) {
								npops2 <- 2 # Moneenko osaan jaetaan
							}
							T2 <- cluster_own(Z2, npops2)
							muutokset <- laskeMuutokset3(
								T2, inds2, rows, data, adjprior, priorTerm, pop
							)
							isoin <- max_MATLAB(muutokset)[[1]]
							indeksi <- max_MATLAB(muutokset)[[2]]
							if (isoin > maxMuutos) {
								maxMuutos <- isoin
								muuttuvaPop2 <- indeksi %% npops2
								if (muuttuvaPop2==0) muuttuvaPop2 <- npops2
								muuttuvat <- inds2[find(T2 == muuttuvaPop2)]
								i2 <- ceiling(indeksi / npops2)
							}
						}
					}
					if (maxMuutos > 1e-5) {
						muutoksia <- 1
						kokeiltu <- zeros(nRoundTypes, 1)
						rivit <- list()
						for (i in 1:length(muuttuvat)) {
							ind <- muuttuvat[i]
							lisa <- rows[ind, 1]:rows[ind, 2]
							rivit <- rbind(rivit, t(lisa))
						}
						diffInCounts <- computeDiffInCounts(
							t(rivit), size(COUNTS, 1), size(COUNTS, 2), data
						)
						i1 <- PARTITION(muuttuvat[1])
						updateGlobalVariables3(
							muuttuvat, diffInCounts, adjprior, priorTerm, i2
						)
						logml <- logml + maxMuutos
						if (dispText) {
							if (round == 3) {
								print('Action 3')
							} else {
								print('Action 4')
							}
						}
						if (logml > worstLogml) {
							partitionSummary_added <- addToSummary(
								logml, partitionSummary, worstIndex
							)
							partitionSummary <- partitionSummary_added$partitionSummary
							added <- partitionSummary_added$added
							if (added==1) {
								worstLogml <- min_MATLAB(partitionSummary[, 2])[[1]]
								worstIndex <- min_MATLAB(partitionSummary[, 2])[[2]]
							}
						}

					} else {
						kokeiltu[round] <- 1
					}
				} else if (round == 5 || round == 6) {
					j <- 0
					muutettu <- 0
					poplogml <- POP_LOGML
					partition <- PARTITION
					counts <- COUNTS
					sumcounts <- SUMCOUNTS
					logdiff <- LOGDIFF

					pops <- sample(npops)
					while (j < npops & muutettu == 0) {
						j <- j + 1
						pop <- pops[j]
						totalMuutos <- 0
						inds <- find(PARTITION == pop)
						if (round == 5) {
							aputaulu <- c(inds, rand(length(inds), 1))
							aputaulu <- sortrows(aputaulu, 2)
							inds <- t(aputaulu[, 1])
						} else if (round == 6) {
							inds <- returnInOrder(
								inds, pop, rows, data, adjprior, priorTerm
							)
						}

						i <- 0

						while (length(inds) > 0 & i < length(inds)) {
							i <- i + 1
							ind <- inds[i]

							muutokset_diffInCounts <- laskeMuutokset(
								ind, rows, data, adjprior, priorTerm
							)
							muutokset <- muutokset_diffInCounts$muutokset
							diffInCounts <- muutokset_diffInCounts$diffInCounts

							muutokset[pop] <- -1e50 # Varmasti ei suurin!!!
							maxMuutos <- max_MATLAB(muutokset)[[1]]
							i2 <- max_MATLAB(muutokset)[[2]]
							updateGlobalVariables(
								ind, i2, diffInCounts, adjprior, priorTerm
							)

							totalMuutos <- totalMuutos+maxMuutos
							logml <- logml + maxMuutos
							if (round == 6) {
								# Lopetetaan heti kun muutos on positiivinen.
								if (totalMuutos > 1e-5) {
									i <- length(inds)
								}
							}
						}

						if (totalMuutos > 1e-5) {
							kokeiltu <- zeros(nRoundTypes, 1)
							muutettu <- 1
							if (muutoksia == 0) {
								muutoksia <- 1  # Ulompi kirjanpito.
								if (dispText) {
									if (round == 5) {
										print('Action 5')
									} else {
										print('Action 6')
									}
								}
							}
							if (logml > worstLogml) {
								partitionSummary_added <- addToSummary(
									logml, partitionSummary, worstIndex
								)
								partitionSummary <- partitionSummary_added$partitionSummary
								added <- partitionSummary_added$added
								if (added==1) {
									worstLogml <- min_MATLAB(partitionSummary[, 2])[[1]]
									worstIndex <- min_MATLAB(partitionSummary[, 2])[[2]]
								}
							}
						} else {
							# Miss��n vaiheessa tila ei parantunut.
							# Perutaan kaikki muutokset.
							PARTITION <- partition
							SUMCOUNTS <- sumcounts
							POP_LOGML <- poplogml
							COUNTS <- counts
							logml <- logml - totalMuutos
							LOGDIFF <- logdiff
							kokeiltu[round] <- 1
						}
					}
					rm(partition, sumcounts, counts, poplogml)

				} else if (round == 7) {
					emptyPop <- findEmptyPop(npops)
					j <- 0
					pops <- sample(npops)
					muutoksiaNyt <- 0
					if (emptyPop == -1) {
						j <- npops
					}
					while (j < npops) {
						j <- j + 1
						pop <- pops[j]
						inds2 <- find(PARTITION == pop)
						ninds2 <- length(inds2)
						if (ninds2 > 5) {
							partition <- PARTITION
							sumcounts <- SUMCOUNTS
							counts <- COUNTS
							poplogml <- POP_LOGML
							logdiff <- LOGDIFF

							dist2 <- laskeOsaDist(inds2, dist, ninds);
							Z2 <- linkage(t(dist2))
							T2 <- cluster_own(Z2, 2)
							muuttuvat <- inds2[find(T2 == 1)]

							muutokset <- laskeMuutokset3(
								T2, inds2, rows, data, adjprior, priorTerm, pop
							)
							totalMuutos <- muutokset(1, emptyPop)

							rivit <- list()
							for (i in 1:length(muuttuvat)) {
								ind <- muuttuvat[i]
								lisa <- rows[ind, 1]:rows[ind, 2]
								rivit <- c(rivit, lisa)
							}
							diffInCounts <- computeDiffInCounts(
								rivit, size(COUNTS, 1), size(COUNTS, 2), data
							)

							updateGlobalVariables3(
								muuttuvat, diffInCounts, adjprior, priorTerm,
								emptyPop
							)

							muutettu <- 1
							while (muutettu == 1) {
								muutettu <- 0
								# Siirret��n yksil�it� populaatioiden v�lill�
								muutokset <- laskeMuutokset5(
									inds2, rows, data, adjprior, priorTerm,
									pop, emptyPop
								)

								maxMuutos <- indeksi <- max_MATLAB(muutokset)

								muuttuva <- inds2(indeksi)
								if (PARTITION(muuttuva) == pop) {
									i2 <- emptyPop
								} else {
									i2 <- pop
								}

								if (maxMuutos > 1e-5) {
									rivit <- rows[muuttuva, 1]:rows[muuttuva, 2]
									diffInCounts <- computeDiffInCounts(
										rivit, size(COUNTS, 1), size(COUNTS, 2),
										data
									)
									updateGlobalVariables3(
										muuttuva,diffInCounts, adjprior,
										priorTerm, i2
									)
									muutettu <- 1
									totalMuutos <- totalMuutos + maxMuutos
								}
							}
							if (totalMuutos > 1e-5) {
								muutoksia <- 1
								logml <- logml + totalMuutos
								if (logml > worstLogml) {
									partitionSummary_added = addToSummary(
										logml, partitionSummary, worstIndex
									)
									partitionSummary_added <- partitionSummary_added$partitionSummary
									added <- partitionSummary_added$added
									if (added == 1) {
										worstLogml <- min_MATLAB(partitionSummary[, 2])[[1]]
										worstIndex <- min_MATLAB(partitionSummary[, 2])[[2]]
									}
								}
								if (muutoksiaNyt == 0) {
									if (dispText) {
										print('Action 7')
									}
									muutoksiaNyt <- 1
								}
								kokeiltu <- zeros(nRoundTypes, 1)
								j <- npops
							} else {
								# palutetaan vanhat arvot
								PARTITION <- partition
								SUMCOUNTS <- sumcounts
								COUNTS <- counts
								POP_LOGML <- poplogml
								LOGDIFF <- logdiff
							}
						}

					}

					if (muutoksiaNyt == 0) {
						kokeiltu[round] <- 1
					}

				}
			}

			if (muutoksia == 0) {
				if (vaihe <= 4) {
					vaihe <= vaihe + 1
				} else if (vaihe == 5) {
					ready <- 1
				}
			} else {
				muutoksia <- 0
			}

			if (ready == 0) {
				if (vaihe == 1) {
					roundTypes <- c(1)
				} else if (vaihe == 2) {
					roundTypes <- c(2, 1)
				} else if (vaihe == 3) {
					roundTypes <- c(5, 5, 7)
				} else if (vaihe == 4) {
					roundTypes = c(4, 3, 1)
				} else if (vaihe == 5) {
					roundTypes <- c(6, 7, 2, 3, 4, 1)
				}
			}
		}

		# TALLENNETAAN

		npops <- poistaTyhjatPopulaatiot(npops)
		POP_LOGML <- computePopulationLogml(1:npops, adjprior, priorTerm)
		if (dispText) {
			print(paste('Found partition with', num2str(npops), 'populations.'))
			print(paste('Log(ml) =', as.character(logml)))
			print(' ')
		}

		if (logml > logmlBest) {
			# P�ivitet��n parasta l�ydetty� partitiota.
			logmlBest <- logml
			npopsBest <- npops
			partitionBest <- PARTITION
			countsBest <- COUNTS
			sumCountsBest <- SUMCOUNTS
			pop_logmlBest <- POP_LOGML
			logdiffbest <- LOGDIFF
		}
	}
	return(
		list(logml = logml, npops = npops, partitionSummary = partitionSummary)
	)
}

# %-----------------------------------------------------------------------


# function [sumcounts, counts, logml] = ...
#     initialPopCounts(data, npops, rows, noalle, adjprior)

# nloci=size(data,2);
# counts = zeros(max(noalle),nloci,npops);
# sumcounts = zeros(npops,nloci);

# for i=1:npops
#     for j=1:nloci
#         i_rivit = rows(i,1):rows(i,2);
#         havainnotLokuksessa = find(data(i_rivit,j)>=0);
#         sumcounts(i,j) = length(havainnotLokuksessa);
#         for k=1:noalle(j)
#             alleleCode = k;
#             N_ijk = length(find(data(i_rivit,j)==alleleCode));
#             counts(k,j,i) = N_ijk;
#         end
#     end
# end

# logml = laskeLoggis(counts, sumcounts, adjprior);

# %-----------------------------------------------------------------------


# function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
# % Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
# % lukum��r�t (vastaavasti kuin COUNTS:issa), jotka ovat data:n
# % riveill� rows. rows pit�� olla vaakavektori.

# diffInCounts = zeros(max_noalle, nloci);
# for i=rows
#     row = data(i,:);
#     notEmpty = find(row>=0);

#     if length(notEmpty)>0
#         diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
#             diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
#     end
# end

# %----------------------------------------------------------------------------


# function dist2 = laskeOsaDist(inds2, dist, ninds)
# % Muodostaa dist vektorista osavektorin, joka sis�lt�� yksil�iden inds2
# % v�liset et�isyydet. ninds=kaikkien yksil�iden lukum��r�.

# ninds2 = length(inds2);
# apu = zeros(nchoosek(ninds2,2),2);
# rivi = 1;
# for i=1:ninds2-1
#     for j=i+1:ninds2
#         apu(rivi, 1) = inds2(i);
#         apu(rivi, 2) = inds2(j);
#         rivi = rivi+1;
#     end
# end
# apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
# dist2 = dist(apu);

# %---------------------------------------------------------

# function T = cluster_own(Z,nclust)
# true=logical(1);
# false=logical(0);
# maxclust = nclust;
# % Start of algorithm
# m = size(Z,1)+1;
# T = zeros(m,1);
#    % maximum number of clusters based on inconsistency
#    if m <= maxclust
#       T = (1:m)';
#    elseif maxclust==1
#       T = ones(m,1);
#    else
#       clsnum = 1;
#       for k = (m-maxclust+1):(m-1)
#          i = Z(k,1); % left tree
#          if i <= m % original node, no leafs
#             T(i) = clsnum;
#             clsnum = clsnum + 1;
#          elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
#             T = clusternum(Z, T, i-m, clsnum);
#             clsnum = clsnum + 1;
#          end
#          i = Z(k,2); % right tree
#          if i <= m  % original node, no leafs
#             T(i) = clsnum;
#             clsnum = clsnum + 1;
#          elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
#             T = clusternum(Z, T, i-m, clsnum);
#             clsnum = clsnum + 1;
#          end
#       end
#    end

# function T = clusternum(X, T, k, c)
# m = size(X,1)+1;
# while(~isempty(k))
#    % Get the children of nodes at this level
#    children = X(k,1:2);
#    children = children(:);

#    % Assign this node number to leaf children
#    t = (children<=m);
#    T(children(t)) = c;

#    % Move to next level
#    k = children(~t) - m;
# end
