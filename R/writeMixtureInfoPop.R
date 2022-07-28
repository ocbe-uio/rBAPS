#' @title Write Mixture Info Pop
#' @description Writes information about the pop mixture
#' @param logml logml
#' @param rows rows
#' @param data data
#' @param adjprior adjprior
#' @param priorTerm priorTerm
#' @param outPutFile outPutFile
#' @param inputFile inputFile
#' @param partitionSummary partitionSummary
#' @param popnames popnames
#' @param fixedK fixedK
#' @export
writeMixtureInfoPop <- function(logml, rows, data, adjprior, priorTerm,
                                outPutFile, inputFile, partitionSummary,
                                popnames, fixedK) {
  ninds <- size(rows, 1)
  npops <- size(COUNTS, 3)
  names <- size(popnames, 1) == ninds # Tarkistetaan ett?nimet viittaavat yksilöihin
  changesInLogml <- vector()
  if (!missing(outPutFile)) {
    fid <- vector()
  }
  cat("RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:\n")
  cat("Data file:", inputFile, "\n")
  cat("Number of clustered groups:", ownNum2Str(ninds), "\n")
  cat("Number of clusters in optimal partition:", ownNum2Str(npops), "\n")
  cat("Log(marginal likelihood) of optimal partition:", ownNum2Str(logml), "\n")
  if (exists("fid")) {
    append(fid, "RESULTS OF GROUP LEVEL MIXTURE ANALYSIS:\n")
    append(fid, c("Data file:", inputFile, "\n"))
    append(fid, c("Number of clustered groups:", ownNum2Str(ninds), "\n"))
    append(fid, c("Number of clusters in optimal partition:", ownNum2Str(npops), "\n"))
    append(fid, c("Log(marginal likelihood) of optimal partition:", ownNum2Str(logml), "\n\n"))
  }
  cluster_count <- length(unique(PARTITION))
  cat("Best Partition:\n")
  if (exists("fid")) {
    append(fid, c("Best partition:\n"))
  }
  for (m in 1:cluster_count) {
    indsInM <- find(PARTITION == m)
    length_of_beginning <- 11 + floor(log10(m))
    cluster_size <- length(indsInM)
    if (names) {
      text <- c("Cluster ", as.character(m), ": {", as.character(popnames[indsInM[1]]))
      for (k in 2:cluster_size) {
        text <- c(text, ", ", as.character(popnames[[indsInM[k]]]))
      }
    } else {
      text <- c("Cluster ", as.character(m), ": {", as.character(indsInM[1]))
      for (k in 2:cluster_size) {
        text <- c(text, ", ", as.character(indsInM[k]))
      }
    }
    text <- c(text, "}")
    while (length(text) > 58) {
      # Take one line and display it.
      new_line <- takeLine(text, 58)
      text <- text[(length(new_line) + 1):length(text)]
      cat(new_line, "\n")
      if (exists("fid")) {
        append(fid, c(new_line, "\n"))
      }
      if (length(text) > 0) {
        text <- c(blanks(length_of_beginning), text)
      } else {
        text <- vector()
      }
    }
    if (!is.null(text)) {
      cat(text, "\n")
      if (exists(fid)) {
        append(fid, c("\n", text, "\n"))
      }
    }
  }
  if (npops > 1) {
    cat("\n\nChanges in log(marginal likelihood) if (group i is moved to cluster j:")
    if (exists("fid")) {
      append(fid, " \n \n", )
      append(fid, "Changes in log(marginal likelihood) if (group i is moved to cluster j:")
    }
    if (names) {
      nameSizes <- zeros(ninds, 1)
      for (i in 1:ninds) {
        nimi <- as.character(popnames[i])
        nameSizes[i] <- length(nimi)
      }
      maxSize <- max(nameSizes)
      maxSize <- max(maxSize, 5)
      erotus <- maxSize - 5
      alku <- blanks(erotus)
      ekarivi <- c(alku, "group", blanks(6 + erotus))
    } else {
      ekarivi <- "group      "
    }
    for (i in 1:cluster_count) {
      ekarivi <- c(ekarivi, ownNum2Str(i), blanks(8 - floor(log10(i))))
    }
    cat(ekarivi, "\n")
    if (exists("fid")) {
      append(fid, c(ekarivi, "\n"))
    }
    changesInLogml <- t(LOGDIFF)
    for (ind in 1:ninds) {
      muutokset <- changesInLogml[, ind]
      if (names) {
        nimi <- as.character(popnames[ind])
        rivi <- c(blanks(maxSize - length(nimi)), nimi, ":")
      } else {
        rivi <- c(blanks(4 - floor(log10(ind))), ownNum2Str(ind), ":")
      }
      for (j in 1:npops) {
        rivi <- c(rivi, "  ", logml2String(omaRound(muutokset(j))))
      }
      cat(rivi, "\n")
      if (exists("fid")) {
        append(fid, c(rivi, "\n"))
      }
    }
    cat("  ")
    cat("KL-divergence matrix in PHYLIP format:")
    dist_mat <- zeros(npops, npops)
    if (exists("fid")) {
      append(fid, " \n")
      append(fid, " \n")
      append(fid, "KL - divergence matrix in PHYLIP format:\n")
    }
    maxnoalle <- size(COUNTS, 1)
    nloci <- size(COUNTS, 2)
    d <- zeros(maxnoalle, nloci, npops)
    prior <- adjprior
    prior[find[prior == 1]] <- 0
    nollia <- find(all(prior == 0)) # Lokukset, joissa oli havaittu vain yht?alleelia.
    prior[1, nollia] <- 1
    for (pop1 in 1:npops) {
      d[, , pop1] <- (squeeze(COUNTS[, , pop1]) + prior) /
        repmat(sum(squeeze(COUNTS[, , pop1]) + prior), c(maxnoalle, 1))
    }
    ekarivi <- as.character(npops)
    cat(ekarivi, "\n")
    if (exists("fid")) {
      append(fid, c(ekarivi, "\n"))
    }
    for (pop1 in 1:npops) {
      rivi <- c(blanks(2 - floor(log10(pop1))), as.character(pop1), "  ")
      for (pop2 in 1:(pop1 - 1)) {
        dist1 <- d[, , pop1]
        dist2 <- d[, , pop2]
        div12 <- sum(sum(dist1 * log2((dist1 + 10^-10) / (dist2 + 10^-10)))) / nloci
        div21 <- sum(sum(dist2 * log2((dist2 + 10^-10) / (dist1 + 10^-10)))) / nloci
        div <- (div12 + div21) / 2
        dist_mat[pop1, pop2] <- div
      }
    }
    dist_mat <- dist_mat + t(dist_mat) # make it symmetric
    for (pop1 in 1:npops) {
      rivi <- c("Cluster_", as.character(pop1), " ")
      for (pop2 in 1:npops) {
        rivi <- c(rivi, kldiv2str(dist_mat(pop1, pop2)), " ")
      }
      cat(rivi)
      if (exists("fid")) {
        append(fid, c(rivi, "\n"))
      }
    }
  }
  cat(" \n \n \n")
  cat(
    "List of sizes of 10 best visited partitions and corresponding",
    "log(ml) values\n"
  )
  if (exists("fid")) {
    append(fid, " \n\n")
    append(fid, " \n\n")
    append(fid, " \n\n")
    append(fid, " \n\n")
    append(
      fid,
      cat(
        "List of sizes of 10 best visited partitions and corresponding",
        "log(ml) values\n"
      )
    )

  }
  partitionSummary <- sortrows(partitionSummary, 2)
  partitionSummary <- partitionSummary[size(partitionSummary, 1):-1, ]
  partitionSummary <- partitionSummary[find(partitionSummary[, 2] > -1e49), ]
  if (size(partitionSummary, 1) > 10) {
    vikaPartitio <- 10
  } else {
    vikaPartitio <- size(partitionSummary, 1)
  }
  for (part in 1:vikaPartitio) {
    line <- c(
      as.character(partitionSummary[part, 1]),
      "    ",
      as.character(partitionSummary[part, 2])
    )
    cat(line, "\n")
    if (exists("fid")) {
      append(fid, c(line, "\n"))
    }
  }
  if (!fixedK) {
    cat(" \n")
    cat(" \n")
    cat("Probabilities for number of clusters\n")
    if (exists("fid")) {
      append(fid, " \n\n")
      append(fid, " \n\n")
      append(fid, "Probabilities for number of clusters\n")
    }
    npopsTaulu <- unique(partitionSummary[, 1])
    len <- length(npopsTaulu)
    probs <- zeros(len, 1)
    partitionSummary[, 2] <- partitionSummary[, 2] - max(partitionSummary[, 2])
    sumtn <- sum(exp(partitionSummary[, 2]))
    for (i in 1:len) {
      npopstn <- sum(
        exp(partitionSummary(find(partitionSummary[, 1] == npopsTaulu(i)), 2))
      )
      probs[i] <- npopstn / sumtn
    }
    for (i in 1:len) {
      if (probs(i) > 1e-5) {
        line <- c(as.character(npopsTaulu(i)), "   ", as.character(probs(i)))
        cat(line)
        if (exists("fid")) {
          append(fid, line)
          append(fid, "\n")
        }
      }
    }
  }
  if (exists("fid")) {
    save(fid, file = outPutFile)
  }
  return(changesInLogml)
}
