#' @title Write Mixture Info
#' @description Writes information about the mixture
#' @param logml logml
#' @param rowsFromInd rowsFromInd
#' @param data data
#' @param adjprior adjprior
#' @param priorTerm priorTerm
#' @param outPutFile outPutFile
#' @param inputFile inputFile
#' @param partitionSummary partitionSummary
#' @param popnames popnames
#' @param fixedK fixedK
#' @export
writeMixtureInfo <- function(logml, rowsFromInd, data, adjprior, priorTerm, outPutFile, inputFile, partitionSummary, popnames, fixedK) {
  changesInLogml <- list()
  ninds <- size(data, 1) / rowsFromInd
  npops <- size(COUNTS, 3)
  # Check that the names refer to individuals
  names <- (size(popnames, 1) == ninds) # Tarkistetaan ett?nimet viittaavat yksil�ihin

  if (length(outPutFile) > 0) {
    fid <- load(outPutFile)
  } else {
    fid <- -1
    # TODO: replace sink with option that will record input and output
    sink("baps4_output.baps", split = TRUE) # save in text anyway.
  }

  dispLine()
  cat("RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:")
  cat(c("Data file: ", inputFile))
  cat("Model: independent")
  cat(c("Number of clustered individuals: ", ownNum2Str(ninds)))
  cat(c("Number of groups in optimal partition: ", ownNum2Str(npops)))
  cat(c("Log(marginal likelihood) of optimal partition: ", ownNum2Str(logml)))
  cat(" ")
  if (fid != -1) {
    append(fid, "RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:\n")
    append(fid, c("Data file: ", inputFile, "\n"))
    append(
      fid,
      c("Number of clustered individuals: ", ownNum2Str(ninds), "\n")
    )
    append(
      fid,
      c(
        "Number of groups in optimal partition: ",
        ownNum2Str(npops), "\n"
      )
    )
    append(
      fid,
      c(
        "Log(marginal likelihood) of optimal partition: ",
        ownNum2Str(logml),
        "\n"
      )
    )
  }

  cluster_count <- length(unique(PARTITION))
  cat("Best Partition: ")
  if (fid != -1) {
    append(fid, c("Best Partition: ", "\n"))
  }
  for (m in 1:cluster_count) {
    indsInM <- find(PARTITION == m)
    length_of_beginning <- 11 + floor(log10(m))
    cluster_size <- length(indsInM)

    if (names) {
      text <- c(
        "Cluster ",
        as.character(m),
        ": {",
        as.character(popnames[[indsInM[1]]])
      )
      for (k in 2:cluster_size) {
        text <- c(text, ", ", as.character(popnames[[indsInM[k]]]))
      }
    } else {
      text <- c(
        "Cluster ", as.character(m), ": {", as.character(indsInM[1])
      )
      for (k in 2:cluster_size) {
        text <- c(text, ", ", as.character(indsInM[k]))
      }
    }
    text <- c(text, "}")
    while (length(text) > 58) {
      # Take one line and display it.
      new_line <- takeLine(text, 58)
      text <- (length(new_line) + 1):length(text)
      cat(new_line)
      if (fid != -1) {
        append(fid, new_line)
        append(fid, "\n")
      }
      if (length(text) > 0) {
        text <- c(blanks(length_of_beginning), text)
      } else {
        text <- ""
      }
    }
    if (text != "") {
      cat(text)
      if (fid != -1) {
        append(fid, text)
        append(fid, "\n")
      }
    }
  }

  if (npops > 1) {
    cat(" ")
    cat(" ")
    cat(
      "Changes in log(marginal likelihood)",
      " if indvidual i is moved to group j:"
    )
    if (fid != -1) {
      append(fid, " ")
      append(fid, "\n")
      append(fid, " ")
      append(fid, "\n")
      append(
        fid,
        c(
          "Changes in log(marginal likelihood)",
          "if indvidual i is moved to group j:"
        )
      )
      append(fid, "\n")
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
      ekarivi <- c(alku, "  ind", blanks(6 + erotus))
    } else {
      ekarivi <- "  ind      "
    }
    for (i in 1:cluster_count) {
      ekarivi <- c(ekarivi, ownNum2Str(i), blanks(8 - floor(log10(i))))
    }
    cat(ekarivi)
    if (fid != -1) {
      append(fid, ekarivi)
      append(fid, "\n")
    }

    # %ninds = size(data,1)/rowsFromInd;
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
        rivi <- c(rivi, "  ", logml2String(omaRound(muutokset[j])))
      }
      cat(rivi)
      if (fid != -1) {
        append(fid, rivi)
        append(fid, "\n")
      }
    }

    cat(" ")
    cat(" ")
    cat("KL-divergence matrix in PHYLIP format:")

    dist_mat <- zeros(npops, npops)
    if (fid != -1) {
      append(fid, " ")
      append(fid, " ")
      append(fid, c("KL-divergence matrix in PHYLIP format:"))
      append(fid, "\n")
    }

    maxnoalle <- size(COUNTS, 1)
    nloci <- size(COUNTS, 2)
    d <- zeros(maxnoalle, nloci, npops)
    prior <- adjprior
    prior[find(prior == 1)] <- 0
    nollia <- find(all(prior == 0)) # Loci in which only one allele was detected.
    prior[1, nollia] <- 1
    for (pop1 in 1:npops) {
      d[, , pop1] <- (squeeze(COUNTS[, , pop1]) + prior) /
        repmat(sum(squeeze(COUNTS[, , pop1]) + prior), c(maxnoalle, 1))
    }
    ekarivi <- as.character(npops)
    cat(ekarivi)
    if (fid != -1) {
      append(fid, ekarivi)
      append(fid, "\n")
    }

    for (pop1 in 1:npops) {
      for (pop2 in 1:(pop1 - 1)) {
        dist1 <- d[, , pop1]
        dist2 <- d[, , pop2]
        div12 <- sum(
          sum(dist1 * log2((dist1 + 10^-10) / (dist2 + 10^-10)))
        ) / nloci
        div21 <- sum(
          sum(dist2 * log2((dist2 + 10^-10) / (dist1 + 10^-10)))
        ) / nloci
        div <- (div12 + div21) / 2
        dist_mat[pop1, pop2] <- div
      }
    }


    dist_mat <- dist_mat + t(dist_mat) # make it symmetric
    for (pop1 in 1:npops) {
      rivi <- c("Cluster_", as.character(pop1), " ")
      for (pop2 in 1:npops) {
        rivi <- c(rivi, kldiv2str(dist_mat[pop1, pop2]), " ")
      }
      cat(rivi)
      if (fid != -1) {
        append(fid, rivi)
        append(fid, "\n")
      }
    }
  }

  cat(" ")
  cat(" ")
  cat(
    "List of sizes of 10 best visited partitions",
    "and corresponding log(ml) values"
  )

  if (fid != -1) {
    append(fid, " ")
    append(fid, "\n")
    append(fid, " ")
    append(fid, "\n")
    append(
      fid,
      c(
        "List of sizes of 10 best visited partitions",
        "and corresponding log(ml) values"
      )
    )
    append(fid, "\n")
  }

  partitionSummary <- sortrows(partitionSummary, 2)
  partitionSummary <- partitionSummary[size(partitionSummary, 1):1, ]
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
      as.character(partitionSummary(part, 2))
    )
    cat(line)
    if (fid != -1) {
      append(fid, line)
      append(fid, "\n")
    }
  }

  if (!fixedK) {
    cat(" ")
    cat(" ")
    cat("Probabilities for number of clusters")

    if (fid != -1) {
      append(fid, " ")
      append(fid, "\n")
      append(fid, " ")
      append(fid, "\n")
      append(fid, c("Probabilities for number of clusters"))
      append(fid, "\n")
    }

    npopsTaulu <- unique(partitionSummary[, 1])
    len <- length(npopsTaulu)
    probs <- zeros(len, 1)
    partitionSummary[, 2] <- partitionSummary[, 2] -
      max(partitionSummary[, 2])
    sumtn <- sum(exp(partitionSummary[, 2]))
    for (i in 1:len) {
      npopstn <- sum(
        exp(
          partitionSummary[find(
            partitionSummary[, 1] == npopsTaulu[i]
          ), 2]
        )
      )
      probs[i] <- npopstn / sumtn
    }
    for (i in 1:len) {
      if (probs[i] > 1e-5) {
        line <- c(
          as.character(npopsTaulu[i]), "   ", as.character(probs[i])
        )
        cat(line)
        if (fid != -1) {
          append(fid, line)
          append(fid, "\n")
        }
      }
    }
  }
  return(changesInLogml)
}
