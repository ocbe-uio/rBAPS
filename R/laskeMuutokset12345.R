#' @title Calculate changes (spatial mixture class)
spatialMixture_muutokset <- R6Class(
  classname = "spatialMixture_muutokset",
  public = list(
    #' @param ind ind
    #' @param rowsFromInd rowsFromInd
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param logml logml
    #' @param cliques cliques
    #' @param separators separators
    laskeMuutokset = function(
      ind, rowsFromInd, data, adjprior, priorTerm, logml, cliques, separators
    ) {
      #  Palauttaa npops * 1 taulun, jossa i:s alkio kertoo, mik?olisi
      #  muutos logml:ss? mikהli yksil?ind siirretההn koriin i.
      #  diffInCounts on poistettava COUNTS:in siivusta i1 ja lisהttהv?
      #  COUNTS:in siivuun i2, mikהli muutos toteutetaan.
      npops <- size(COUNTS, 3)
      muutokset <- zeros(npops, 1)

      emptyPop_pops <- findEmptyPop(npops)
      emptyPop <- emptyPop_pops$emptyPop
      pops <- emptyPop_pops$pops
      rm(emptyPop_pops)

      i1 <- PARTITION(ind)
      i2 <- pops[find(pops != i1)]
      if (emptyPop > 0) {
        i2 <- c(i2, emptyPop)
      }

      rows <- ((ind - 1) * rowsFromInd + 1):(ind * rowsFromInd)
      diffInCounts <- computeDiffInCounts(rows, size(COUNTS, 1), size(COUNTS, 2), data)
      diffInSumCounts <- sum(diffInCounts)

      diffInCliqCounts <- computeDiffInCliqCounts(cliques, ind)
      diffInSepCounts <- computeDiffInCliqCounts(separators, ind)

      COUNTS[, ,i1] <- COUNTS[, , i1] - diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
      CLIQCOUNTS[, i1] <- CLIQCOUNTS[, i1] - diffInCliqCounts
      SEPCOUNTS[, i1] <- SEPCOUNTS[, i1] - diffInSepCounts

      for (i in i2) {
        CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] + diffInCliqCounts
        SEPCOUNTS[, i] <- SEPCOUNTS[, i] + diffInSepCounts
        COUNTS[, ,i] <- COUNTS[, , i] + diffInCounts
        SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] + diffInSumCounts

        muutokset[i] <- computeLogml(adjprior, priorTerm) - logml

        CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] - diffInCliqCounts
        SEPCOUNTS[, i] <- SEPCOUNTS[, i] - diffInSepCounts
        COUNTS[, , i] <- COUNTS[, , i] - diffInCounts
        SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] - diffInSumCounts
      }

      COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts
      CLIQCOUNTS[, i1] <- CLIQCOUNTS[, i1] + diffInCliqCounts
      SEPCOUNTS[, i1] <- SEPCOUNTS[, i1] + diffInSepCounts

      #  Asetetaan muillekin tyhjille populaatioille sama muutos, kuin
      #  emptyPop:lle

      if (emptyPop > 0) {
        empties <- mysetdiff(1:npops, c(i2, i1))
        muutokset[empties] <- muutokset(emptyPop)
      }
      return(list(muutokset = muutokset, diffInCounts = diffInCounts))
    },
    #' @param i1 i1
    #' @param rowsFromInd rowsFromInd
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param logml logml
    #' @param cliques cliques
    #' @param separators separators
    laskeMuutokset2 = function(
      i1, rowsFromInd, data, adjprior, priorTerm, logml, cliques, separators
    ) {
      #  Palauttaa npops * 1 taulun, jossa i:s alkio kertoo, mik?olisi
      #  muutos logml:ss? mikהli korin i1 kaikki yksilצt siirretההn
      #  koriin i.
      #  Laskee muutokset vain yhdelle tyhjהlle populaatiolle, muille tulee
      #  muutokseksi 0.
      # global COUNTS      # global SUMCOUNTS
      # global PARTITION   # global POP_LOGML
      # global CLIQCOUNTS  # global SEPCOUNTS

      npops <- size(COUNTS, 3)
      muutokset <- zeros(npops, 1)

      emptyPop <- findEmptyPop(npops)$emptyPop
      pops <- findEmptyPop(npops)$npops

      i2 <- pops[find(pops != i1)]
      if (emptyPop > 0) {
        i2 <- c(i2, emptyPop)
      }

      inds <- find(PARTITION == i1)
      rows <- computeRows(rowsFromInd, inds, length(inds))

      diffInCounts <- computeDiffInCounts(rows, size(COUNTS, 1), size(COUNTS, 2), data)
      diffInSumCounts <- sum(diffInCounts)
      diffInCliqCounts <- computeDiffInCliqCounts(cliques, inds)
      diffInSepCounts <- computeDiffInCliqCounts(separators, inds)

      COUNTS[, ,i1] <- COUNTS[, , i1] - diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
      CLIQCOUNTS[, i1] <- 0
      SEPCOUNTS[, i1] <- 0

      for (i in i2) {
        CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] + diffInCliqCounts
        SEPCOUNTS[, i] <- SEPCOUNTS[, i] + diffInSepCounts
        COUNTS[, ,i] <- COUNTS[, , i] + diffInCounts
        SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] + diffInSumCounts

        muutokset[i] <- computeLogml(adjprior, priorTerm) - logml

        CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] - diffInCliqCounts
        SEPCOUNTS[, i] <- SEPCOUNTS[, i] - diffInSepCounts
        COUNTS[, ,i] <- COUNTS[, , i] - diffInCounts
        SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] - diffInSumCounts
      }

      COUNTS[, ,i1] <- COUNTS[, , i1] + diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts
      CLIQCOUNTS[, i1] <- diffInCliqCounts
      SEPCOUNTS[, i1] <- diffInSepCounts
      return(list(muutokset = muutokset, diffInCounts = diffInCounts))
    },
    #' @param T2 T2
    #' @param inds2 inds2
    #' @param rowsFromInd rowsFromInd
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param i1 i1
    #' @param logml logml
    #' @param cliques cliques
    #' @param separators separators
    laskeMuutokset3 = function(
      T2, inds2, rowsFromInd, data, adjprior, priorTerm, i1, logml, cliques,
      separators
    ) {
      #  Palauttaa length(unique(T2)) * npops taulun, jossa (i, j):s alkio
      #  kertoo, mik?olisi muutos logml:ss? jos populaation i1 osapopulaatio
      #  inds2(find(T2 == i)) siirretההn koriin j.
      #  Laskee vain yhden tyhjהn populaation, muita kohden muutokseksi jהה 0.

      # global COUNTS      # global SUMCOUNTS
      # global PARTITION   # global POP_LOGML
      # global CLIQCOUNTS  # global SEPCOUNTS

      npops <- size(COUNTS, 3)
      npops2 <- length(unique(T2))
      muutokset <- zeros(npops2, npops)

      for (pop2 in 1:npops2) {
        inds <- inds2[find(T2 == pop2)]
        ninds <- length(inds)
        if (ninds > 0) {
          rows <- computeRows(rowsFromInd, inds, ninds)

          diffInCounts <- computeDiffInCounts(rows, size(COUNTS, 1), size(COUNTS, 2), data)
          diffInSumCounts <- sum(diffInCounts)
          diffInCliqCounts <- computeDiffInCliqCounts(cliques, inds)
          diffInSepCounts <- computeDiffInCliqCounts(separators, inds)

          COUNTS[, ,i1] <- COUNTS[, , i1] - diffInCounts
          SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
          CLIQCOUNTS[, i1] <- CLIQCOUNTS[, i1] - diffInCliqCounts
          SEPCOUNTS[, i1] <- SEPCOUNTS[, i1] - diffInSepCounts

          emptyPop <- findEmptyPop(npops)$emptyPop
          pops <- findEmptyPop(npops)$pops
          i2 <- pops[find(pops != i1)]
          if (emptyPop > 0) {
            i2 <- c(i2, emptyPop)
          }

          for (i in i2) {
            CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] + diffInCliqCounts
            SEPCOUNTS[, i] <- SEPCOUNTS[, i] + diffInSepCounts
            COUNTS[, ,i] <- COUNTS[, , i] + diffInCounts
            SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] + diffInSumCounts

            muutokset[pop2, i] <- computeLogml(adjprior, priorTerm) - logml

            CLIQCOUNTS[, i] <- CLIQCOUNTS[, i] - diffInCliqCounts
            SEPCOUNTS[, i] <- SEPCOUNTS[, i] - diffInSepCounts
            COUNTS[, ,i] <- COUNTS[, , i] - diffInCounts
            SUMCOUNTS[i, ] <- SUMCOUNTS[i, ] - diffInSumCounts
          }

          COUNTS[, ,i1] <- COUNTS[, , i1] + diffInCounts
          SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts
          CLIQCOUNTS[, i1] <- CLIQCOUNTS[, i1] + diffInCliqCounts
          SEPCOUNTS[, i1] <- SEPCOUNTS[, i1] + diffInSepCounts
        }
      }
      return(muutokset)
    },
    #' @param inds inds
    #' @param rowsFromInd rowsFromInd
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param logml logml
    #' @param cliques cliques
    #' @param separators separators
    #' @param i1 i1
    #' @param i2 i2
    laskeMuutokset5 = function(
      inds, rowsFromInd, data, adjprior, priorTerm, logml, cliques, separators,
      i1, i2
    ) {
      #  Palauttaa length(inds) * 1 taulun, jossa i:s alkio kertoo, mik?olisi
      #  muutos logml:ss? mikהli yksil?i vaihtaisi koria i1:n ja i2:n vהlill?

      # global COUNTS    # global SUMCOUNTS
      # global PARTITION
      # global CLIQCOUNTS  # global SEPCOUNTS

      ninds <- length(inds)
      muutokset <- zeros(ninds, 1)
      cliqsize <- size(CLIQCOUNTS, 2)
      sepsize <- size(SEPCOUNTS, 2)

      for (i in 1:ninds) {
        ind <- inds[i]
        if (PARTITION[ind] == i1) {
          pop1 <- i1 # mist?
          pop2 <- i2 # mihin
        } else {
          pop1 <- i2
          pop2 <- i1
        }
        rows <- ((ind - 1) * rowsFromInd + 1):(ind * rowsFromInd)

        diffInCounts <- computeDiffInCounts(rows, size(COUNTS, 1), size(COUNTS, 2), data)
        diffInSumCounts <- sum(diffInCounts)
        diffInCliqCounts <- computeDiffInCliqCounts(cliques, ind)
        diffInSepCounts <- computeDiffInCliqCounts(separators, ind)

        COUNTS[, ,pop1] <- COUNTS[, , pop1] - diffInCounts
        SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] - diffInSumCounts
        COUNTS[, ,pop2] <- COUNTS[, , pop2] + diffInCounts
        SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] + diffInSumCounts

        CLIQCOUNTS[, pop1] <- CLIQCOUNTS[, pop1] - diffInCliqCounts
        CLIQCOUNTS[, pop2] <- CLIQCOUNTS[, pop2] + diffInCliqCounts
        SEPCOUNTS[, pop1] <- SEPCOUNTS[, pop1] - diffInSepCounts
        SEPCOUNTS[, pop2] <- SEPCOUNTS[, pop2] + diffInSepCounts

        muutokset[i] <- computeLogml(adjprior, priorTerm) - logml

        COUNTS[, ,pop1] <- COUNTS[, , pop1] + diffInCounts
        SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] + diffInSumCounts
        COUNTS[, ,pop2] <- COUNTS[, , pop2] - diffInCounts
        SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] - diffInSumCounts

        CLIQCOUNTS[, pop1] <- CLIQCOUNTS[, pop1] + diffInCliqCounts
        CLIQCOUNTS[, pop2] <- CLIQCOUNTS[, pop2] - diffInCliqCounts
        SEPCOUNTS[, pop1] <- SEPCOUNTS[, pop1] + diffInSepCounts
        SEPCOUNTS[, pop2] <- SEPCOUNTS[, pop2] - diffInSepCounts

      }
      return(muutokset)
    }
  )
)

#' @title Calculate changes (admix1 class)
#' @description Palauttaa npops*npops taulun, jonka alkio (i,j) kertoo, mik?on
#' muutos logml:ss? mikäli populaatiosta i siirretään osuuden verran
#' todennäköisyysmassaa populaatioon j. Mikäli populaatiossa i ei ole mitään
#' siirrettävää, on vastaavassa kohdassa rivi nollia.
#' @importFrom R6 R6Class
admix1_muutokset <- R6Class(
  classname = "admix1_muutokset",
  public = list(
    #' @param osuus Percentages?
    #' @param osuusTaulu Percentage table?
    #' @param omaFreqs own Freqs?
    #' @param logml log maximum likelihood
    laskeMuutokset4 = function(osuus, osuusTaulu, omaFreqs, logml) {
      if (isGlobalEmpty(COUNTS)) {
        npops <- 1
      } else {
        npops <- ifelse(is.na(dim(COUNTS)[3]), 1, dim(COUNTS)[3])
      }
      notEmpty <- which(osuusTaulu > 0.005)
      muutokset <- zeros(npops)
      empties <- !notEmpty

      for (i1 in notEmpty) {
        osuusTaulu[i1] <- osuusTaulu[i1] - osuus
        for (i2 in c(colon(1, i1 - 1), colon(i1 + 1, npops))) {
          osuusTaulu[i2] <- osuusTaulu[i2] + osuus
          loggis <- computeIndLogml(omaFreqs, osuusTaulu)

          # Work around Matlab OOB bug
          if (i1 > nrow(muutokset)) {
            muutokset <- rbind(muutokset, muutokset * 0)
          }
          if (i2 > ncol(muutokset)) {
            muutokset <- cbind(muutokset, muutokset * 0)
          }

          muutokset[i1, i2] <- loggis - logml
          osuusTaulu[i2] <- osuusTaulu[i2] - osuus
        }
        osuusTaulu[i1] <- osuusTaulu[i1] + osuus
      }
      return(muutokset)
    }
  )
)

#' @title Calculate changes (greedyMix class)
#' @description Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik� olisi
#' muutos logml:ss�, mik�li yksil� ind siirret��n koriin i.
#' diffInCounts on poistettava COUNTS:in siivusta i1 ja lis�tt�v�
#' COUNTS:in siivuun i2, mik�li muutos toteutetaan.
#'
#' Lis�ys 25.9.2007:
#' Otettu k�ytt��n globaali muuttuja LOGDIFF, johon on tallennettu muutokset
#' logml:ss� siirrett�ess� yksil�it� toisiin populaatioihin.
greedyMix_muutokset <- R6Class(
  classname = "greedyMix_muutokset",
  public = list(
    #' @param ind ind
    #' @param globalRows globalRows
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    laskeMuutokset = function(ind, globalRows, data, adjprior, priorTerm) {
      npops <- size(COUNTS, 3)
      muutokset <- LOGDIFF[ind, ]

      i1 <- PARTITION[ind]
      i1_logml <- POP_LOGML[i1]
      muutokset[i1] <- 0

      if (is.null(dim(globalRows))) {
        rows <- globalRows[1]:globalRows[2]
      } else {
        rows <- globalRows[ind, 1]:globalRows[ind, 2]
      }
      diffInCounts <- computeDiffInCounts(
        rows, size(COUNTS, 1), size(COUNTS, 2), data
      )
      diffInSumCounts <- colSums(diffInCounts)
      COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
      new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
      COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

      i2 <- matlab2r::find(muutokset == -Inf) # Etsit��n populaatiot jotka muuttuneet viime kerran j�lkeen. (Searching for populations that have changed since the last time)
      i2 <- setdiff(i2, i1)
      i2_logml <- POP_LOGML[i2]

      ni2 <- length(i2)

      COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, ni2))
      SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(ni2, 1))
      new_i2_logml <- computePopulationLogml(i2, adjprior, priorTerm)
      COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, ni2))
      SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(ni2, 1))

      muutokset[i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
      LOGDIFF[ind, ] <- muutokset
      return(list(muutokset = muutokset, diffInCounts = diffInCounts))
    },
    #' @param i1 i1
    #' @param globalRows globalRows
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    laskeMuutokset2 = function(i1, globalRows, data, adjprior, priorTerm) {
      # % Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mik� olisi
      # % muutos logml:ss�, mik�li korin i1 kaikki yksil�t siirret��n
      # % koriin i.

      npops <- size(COUNTS, 3)
      muutokset <- zeros(npops, 1)

      i1_logml <- POP_LOGML[i1]

      inds <- matlab2r::find(PARTITION == i1)
      ninds <- length(inds)

      if (ninds == 0) {
        diffInCounts <- zeros(size(COUNTS, 1), size(COUNTS, 2))
        return()
      }

      rows <- list()
      for (i in 1:ninds) {
        ind <- inds(i)
        lisa <- globalRows(ind, 1):globalRows(ind, 2)
        rows <- c(rows, t(lisa))
      }

      diffInCounts <- computeDiffInCounts(
        t(rows), size(COUNTS, 1), size(COUNTS, 2), data
      )
      diffInSumCounts <- sum(diffInCounts)

      COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
      new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
      COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
      SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

      i2 <- c(1:i1 - 1, i1 + 1:npops)
      i2_logml <- POP_LOGML[i2]

      COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, npops - 1))
      SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(npops - 1, 1))
      new_i2_logml <- computePopulationLogml(i2, adjprior, priorTerm)
      COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, npops - 1))
      SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(npops - 1, 1))

      muutokset[i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
      return(list(muutokset = muutokset, diffInCounts = diffInCounts))
    },
    #' @param T2 T2
    #' @param inds2 inds2
    #' @param globalRows globalRows
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param i1 i1
    laskeMuutokset3 = function(
      T2, inds2, globalRows, data, adjprior, priorTerm, i1
    ) {
      # Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
      # kertoo, mik� olisi muutos logml:ss�, jos populaation i1 osapopulaatio
      # inds2(matlab2r::find(T2==i)) siirret��n koriin j.

      npops <- size(COUNTS, 3)
      npops2 <- length(unique(T2))
      muutokset <- zeros(npops2, npops)

      i1_logml <- POP_LOGML[i1]
      for (pop2 in 1:npops2) {
        inds <- inds2[matlab2r::find(T2 == pop2)]
        ninds <- length(inds)
        if (ninds > 0) {
          rows <- list()
          for (i in 1:ninds) {
            ind <- inds[i]
            lisa <- globalRows[ind, 1]:globalRows[ind, 2]
            rows <- c(rows, t(lisa))
          }
          diffInCounts <- computeDiffInCounts(
            t(rows), size(COUNTS, 1), size(COUNTS, 2), data
          )
          diffInSumCounts <- sum(diffInCounts)

          COUNTS[, , i1] <- COUNTS[, , i1] - diffInCounts
          SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] - diffInSumCounts
          new_i1_logml <- computePopulationLogml(i1, adjprior, priorTerm)
          COUNTS[, , i1] <- COUNTS[, , i1] + diffInCounts
          SUMCOUNTS[i1, ] <- SUMCOUNTS[i1, ] + diffInSumCounts

          i2 <- c(1:i1 - 1, i1 + 1:npops)
          i2_logml <- t(POP_LOGML[i2])

          COUNTS[, , i2] <- COUNTS[, , i2] + repmat(diffInCounts, c(1, 1, npops - 1))
          SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] + repmat(diffInSumCounts, c(npops - 1, 1))
          new_i2_logml <- t(computePopulationLogml(i2, adjprior, priorTerm))
          COUNTS[, , i2] <- COUNTS[, , i2] - repmat(diffInCounts, c(1, 1, npops - 1))
          SUMCOUNTS[i2, ] <- SUMCOUNTS[i2, ] - repmat(diffInSumCounts, c(npops - 1, 1))

          muutokset[pop2, i2] <- new_i1_logml - i1_logml + new_i2_logml - i2_logml
        }
      }
      return(muutokset)
    },
    #' @param inds inds
    #' @param globalRows globalRows
    #' @param data data
    #' @param adjprior adjprior
    #' @param priorTerm priorTerm
    #' @param i1 i1
    #' @param i2 i2
    laskeMuutokset5 = function(inds, globalRows, data, adjprior, priorTerm, i1, i2) {
      # Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mik� olisi
      # muutos logml:ss�, mik�li yksil� i vaihtaisi koria i1:n ja i2:n v�lill�.

      ninds <- length(inds)
      muutokset <- zeros(ninds, 1)

      i1_logml <- POP_LOGML[i1]
      i2_logml <- POP_LOGML[i2]

      for (i in 1:ninds) {
        ind <- inds[i]
        if (PARTITION[ind] == i1) {
          pop1 <- i1 # mist�
          pop2 <- i2 # mihin
        } else {
          pop1 <- i2
          pop2 <- i1
        }
        rows <- globalRows[ind, 1]:globalRows[ind, 2]
        diffInCounts <- computeDiffInCounts(
          rows, size(COUNTS, 1), size(COUNTS, 2), data
        )
        diffInSumCounts <- sum(diffInCounts)



        COUNTS[, , pop1] <- COUNTS[, , pop1] - diffInCounts
        SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] - diffInSumCounts
        COUNTS[, , pop2] <- COUNTS[, , pop2] + diffInCounts
        SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] + diffInSumCounts

        new_logmls <- computePopulationLogml(c(i1, i2), adjprior, priorTerm)
        muutokset[i] <- sum(new_logmls)

        COUNTS[, , pop1] <- COUNTS[, , pop1] + diffInCounts
        SUMCOUNTS[pop1, ] <- SUMCOUNTS[pop1, ] + diffInSumCounts
        COUNTS[, , pop2] <- COUNTS[, , pop2] - diffInCounts
        SUMCOUNTS[pop2, ] <- SUMCOUNTS[pop2, ] - diffInSumCounts
      }

      muutokset <- muutokset - i1_logml - i2_logml
      return(muutokset)
    }
  )
)
