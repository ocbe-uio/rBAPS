indMix <- function(c, npops, counts = NULL, sumcounts = NULL, max_iter = 100L, dispText = FALSE) {
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

  if (isfield(c, "dist")) {
    dist <- c$dist
    Z <- c$Z
  } else {
    Z <- NULL
  }


  rm(c)
  nargin <- length(as.list(match.call())) - 1
  if (nargin < 2) {
    dispText <- 1
    npopstext <- matrix()
    ready <- FALSE
    teksti <- "Input upper bound to the number of populations (possibly multiple values)"
    while (!ready) {
      npopstextExtra <- inputdlg(teksti, 1, "20")
      if (isempty(npopstextExtra)) { # Painettu Cancel:ia
        warnings("Empty value provided")
        return()
      }
      npopstextExtra <- npopstextExtra[[1]]
      if (length(npopstextExtra) >= 255) {
        npopstextExtra <- npopstextExtra[1:255]
        npopstext <- c(npopstext, " ", npopstextExtra)
        teksti <- "The input field length limit (255 characters) was reached. Input more values: "
      } else {
        npopstext <- as.numeric(strsplit(as.character(npopstextExtra), " ")[[1]])
        ready <- TRUE
      }
    }
    rm(ready, teksti)
    if (isempty(npopstext) | length(npopstext) == 1) {
      warning("Empty or 1-length vector provided")
      return()
    } else {
      npopsTaulu <- as.numeric(npopstext)
      ykkoset <- matlab2r::find(npopsTaulu == 1)
      npopsTaulu[ykkoset] <- NA # Mik�li ykk�si� annettu yl�rajaksi, ne poistetaan (if ones are given as an upper limit, they are deleted)
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
  data <- data[, 1:(ncol(data) - 1)]

  logmlBest <- -1e50
  partitionSummary <- -1e50 * ones(30, 2) # Tiedot 30 parhaasta partitiosta (npops ja logml)
  partitionSummary[, 1] <- zeros(30, 1)
  worstLogml <- -1e50
  worstIndex <- 1
  for (run in 1:nruns) {
    npops <- npopsTaulu[[run]]
    if (dispText) {
      dispLine()
      cat(
        "Run ", as.character(run), "/", as.character(nruns),
        ", maximum number of populations ", as.character(npops),
        ".\n",
        sep = ""
      )
    }
    ninds <- size(rows, 1)

    if (!is.null(Z)) {
      initialPartition <- admixture_initialization(initData, npops, Z)
      sumcounts_counts_logml <- initialCounts(
        initialPartition, data, npops, rows, noalle, adjprior
      )
      sumcounts <- sumcounts_counts_logml$sumcounts
      counts <- sumcounts_counts_logml$counts
      logml <- sumcounts_counts_logml$logml

      PARTITION <- zeros(ninds, 1)
      for (i in 1:ninds) {
        apu <- rows[i]
        PARTITION[i] <- initialPartition[apu[1]]
      }
    }

    COUNTS <- counts
    SUMCOUNTS <- sumcounts
    POP_LOGML <- computePopulationLogml(seq_len(npops), adjprior, priorTerm)
    LOGDIFF <- repmat(-Inf, c(ninds, npops))

    # PARHAAN MIXTURE-PARTITION ETSIMINEN
    nRoundTypes <- 7
    kokeiltu <- zeros(nRoundTypes, 1)
    roundTypes <- c(1, 1) # Ykk�svaiheen sykli kahteen kertaan.
    ready <- 0
    vaihe <- 1

    if (dispText) {
      message(
        paste0(
          "\nMixture analysis started with initial ",
          as.character(npops), " populations."
        )
      )
    }

    iter <- 1L
    while (ready != 1) {
      iter <- iter + 1L
      if (iter > max_iter) {
        warning("max_iter reached. Stopping.")
        break
      }
      # FIXME: loop caught in here
      muutoksia <- 0

      if (dispText) {
        message(paste("Performing steps:", as.character(roundTypes)))
      }

      for (n in 1:length(roundTypes)) {
        round <- roundTypes[n]
        kivaluku <- 0

        if (kokeiltu[round] == 1) { # Askelta kokeiltu viime muutoksen j�lkeen
        } else if (round == 0 | round == 1) { # Yksil�n siirt�minen toiseen populaatioon.
          inds <- 1:ninds
          aputaulu <- cbind(inds, rand(ninds, 1))
          aputaulu <- sortrows(aputaulu, 2)
          inds <- t(aputaulu[, 1])
          muutosNyt <- 0

          for (ind in inds) {
            i1 <- PARTITION[ind]
            muutokset_diffInCounts <- greedyMix_muutokset$new()
            # FIXME: using 100-length global variables instead of the ones in this function
            muutokset_diffInCounts <- muutokset_diffInCounts$laskeMuutokset(
              ind, rows, data, adjprior, priorTerm
            )
            muutokset <- muutokset_diffInCounts$muutokset
            diffInCounts <- muutokset_diffInCounts$diffInCounts

            if (round == 1) {
              maxMuutos <- matlab2r::max(muutokset)$max
              i2 <- matlab2r::max(muutokset)$idx
            }

            if (i1 != i2 & maxMuutos > 1e-5) {
              # Tapahtui muutos
              muutoksia <- 1
              if (muutosNyt == 0) {
                muutosNyt <- 1
                if (dispText) message("Action 1")
              }
              kokeiltu <- zeros(nRoundTypes, 1)
              kivaluku <- kivaluku + 1
              updateGlobalVariables(
                ind, i2, diffInCounts, adjprior, priorTerm
              )
              logml <- logml + maxMuutos
              if (logml > worstLogml) {
                temp_addToSum <- addToSummary(
                  logml, partitionSummary, worstIndex
                )
                partitionSummary <- temp_addToSum$partitionSummary
                added <- temp_addToSum$added
                if (added == 1) {
                  temp_minMATLAB <- matlab2r::min(
                    partitionSummary[, 2]
                  )
                  worstLogml <- temp_minMATLAB$mins
                  worstIndex <- temp_minMATLAB$idx
                }
              }
            }
          }

          if (muutosNyt == 0) {
            kokeiltu[round] <- 1
          }
        } else if (round == 2) { # Populaation yhdist�minen toiseen.
          maxMuutos <- 0
          for (pop in seq_len(npops)) {
            muutokset_diffInCounts <- greedyMix_muutokset$new()
            muutokset_diffInCounts <- muutokset_diffInCounts$laskeMuutokset2(
              pop, rows, data, adjprior, priorTerm
            )
            muutokset <- muutokset_diffInCounts$muutokset
            diffInCounts <- muutokset_diffInCounts$diffInCounts
            isoin <- matlab2r::max(muutokset)[[1]]
            indeksi <- matlab2r::max(muutokset)[[2]]
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
              cat("Action 2")
            }
            if (logml > worstLogml) {
              temp_addToSum <- addToSummary(
                logml, partitionSummary, worstIndex
              )
              partitionSummary <- temp_addToSum$partitionSummary
              added <- temp_addToSum$added
              if (added == 1) {
                worstLogml <- matlab2r::min(partitionSummary[, 2])[[1]]
                worstIndex <- matlab2r::min(partitionSummary[, 2])[[2]]
              }
            }
          } else {
            kokeiltu[round] <- 1
          }
        } else if (round == 3 || round == 4) { # Populaation jakaminen osiin.
          maxMuutos <- 0
          ninds <- size(rows, 1)
          for (pop in seq_len(npops)) {
            inds2 <- matlab2r::find(PARTITION == pop)
            ninds2 <- length(inds2)
            if (ninds2 > 2) {
              dist2 <- laskeOsaDist(inds2, dist, ninds)
              Z2 <- linkage(t(dist2))
              if (round == 3) {
                npops2 <- base::max(base::min(20, floor(ninds2 / 5)), 2)
              } else if (round == 4) {
                npops2 <- 2 # Moneenko osaan jaetaan
              }
              T2 <- cluster_own(Z2, npops2)
              muutokset_diffInCounts <- greedyMix_muutokset$new()
              muutokset <- muutokset_diffInCounts$laskeMuutokset3(
                T2, inds2, rows, data, adjprior, priorTerm, pop
              )
              isoin <- matlab2r::max(muutokset)[[1]]
              indeksi <- matlab2r::max(muutokset)[[2]]
              if (isoin > maxMuutos) {
                maxMuutos <- isoin
                muuttuvaPop2 <- indeksi %% npops2
                if (muuttuvaPop2 == 0) muuttuvaPop2 <- npops2
                muuttuvat <- inds2[matlab2r::find(T2 == muuttuvaPop2)]
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
                cat("Action 3")
              } else {
                cat("Action 4")
              }
            }
            if (logml > worstLogml) {
              temp_addToSum <- addToSummary(
                logml, partitionSummary, worstIndex
              )
              partitionSummary <- temp_addToSum$partitionSummary
              added <- temp_addToSum$added
              if (added == 1) {
                worstLogml <- matlab2r::min(partitionSummary[, 2])[[1]]
                worstIndex <- matlab2r::min(partitionSummary[, 2])[[2]]
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
            inds <- matlab2r::find(PARTITION == pop)
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
              muutokset_diffInCounts <- greedyMix_muutokset$new()
              muutokset_diffInCounts <- muutokset_diffInCounts$laskeMuutokset(
                ind, rows, data, adjprior, priorTerm
              )
              muutokset <- muutokset_diffInCounts$muutokset
              diffInCounts <- muutokset_diffInCounts$diffInCounts

              muutokset[pop] <- -1e50 # Varmasti ei suurin!!!
              maxMuutos <- matlab2r::max(muutokset)[[1]]
              i2 <- matlab2r::max(muutokset)[[2]]
              updateGlobalVariables(
                ind, i2, diffInCounts, adjprior, priorTerm
              )

              totalMuutos <- totalMuutos + maxMuutos
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
                muutoksia <- 1 # Ulompi kirjanpito.
                if (dispText) {
                  if (round == 5) {
                    cat("Action 5")
                  } else {
                    cat("Action 6")
                  }
                }
              }
              if (logml > worstLogml) {
                temp_addToSum <- addToSummary(
                  logml, partitionSummary, worstIndex
                )
                partitionSummary <- temp_addToSum$partitionSummary
                added <- temp_addToSum$added
                if (added == 1) {
                  worstLogml <- matlab2r::min(partitionSummary[, 2])[[1]]
                  worstIndex <- matlab2r::min(partitionSummary[, 2])[[2]]
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
            inds2 <- matlab2r::find(PARTITION == pop)
            ninds2 <- length(inds2)
            if (ninds2 > 5) {
              partition <- PARTITION
              sumcounts <- SUMCOUNTS
              counts <- COUNTS
              poplogml <- POP_LOGML
              logdiff <- LOGDIFF

              dist2 <- laskeOsaDist(inds2, dist, ninds)
              Z2 <- linkage(t(dist2))
              T2 <- cluster_own(Z2, 2)
              muuttuvat <- inds2[matlab2r::find(T2 == 1)]
              muutokset <- greedyMix_muutokset$new()
              muutokset <- muutokset$laskeMuutokset3(
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
                muutokset <- greedyMix_muutokset$new()
                muutokset <- muutokset$laskeMuutokset5(
                  inds2, rows, data, adjprior, priorTerm,
                  pop, emptyPop
                )

                maxMuutos <- indeksi <- matlab2r::max(muutokset)

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
                    muuttuva, diffInCounts, adjprior,
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
                  temp_addToSum <- addToSummary(
                    logml, partitionSummary, worstIndex
                  )
                  partitionSummary <- temp_addToSum$partitionSummary
                  added <- temp_addToSum$added
                  if (added == 1) {
                    worstLogml <- matlab2r::min(partitionSummary[, 2])[[1]]
                    worstIndex <- matlab2r::min(partitionSummary[, 2])[[2]]
                  }
                }
                if (muutoksiaNyt == 0) {
                  if (dispText) {
                    cat("Action 7")
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
      # FIXME: muutoksia is never 0, so vaihe never equals 5 and ready 1
      print(paste("i1 =", i1, "i2 =", i2, "maxMuutos =", maxMuutos)) # TEMP
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
          roundTypes <- 1
        } else if (vaihe == 2) {
          roundTypes <- c(2, 1)
        } else if (vaihe == 3) {
          roundTypes <- c(5, 5, 7)
        } else if (vaihe == 4) {
          roundTypes <- c(4, 3, 1)
        } else if (vaihe == 5) {
          roundTypes <- c(6, 7, 2, 3, 4, 1)
        }
      }
    }

    # TALLENNETAAN

    npops <- poistaTyhjatPopulaatiot(npops)
    POP_LOGML <- computePopulationLogml(seq_len(npops), adjprior, priorTerm)
    if (dispText) {
      message("Found partition with ", as.character(npops), " populations.")
      message("Log(ml) = ", as.character(logml))
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
