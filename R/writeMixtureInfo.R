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
#' @param verbose if \code{TRUE}, prints extra output information
writeMixtureInfo <- function(
  logml, rowsFromInd, data, adjprior, priorTerm, outPutFile, inputFile,
  partitionSummary, popnames, fixedK, verbose
) {
  ninds <- size(data, 1) / rowsFromInd
  npops <- size(COUNTS, 3)
  # Check that the names refer to individuals

  # Tarkistetaan ett?nimet viittaavat yksil�ihin
  names <- (size(popnames, 1) == ninds)

  if (length(outPutFile) > 0) {
    fid <- load(outPutFile)
  } else {
    fid <- -1
    outPutFile <- file.path(tempdir(), "baps4_output.baps")
    message("Output saved to", outPutFile)
    sink(outPutFile, split = TRUE) # save in text anyway.
  }
  if (verbose) {
    dispLine()
    cat("RESULTS OF INDIVIDUAL LEVEL MIXTURE ANALYSIS:\n")
    cat("Data file: ", inputFile, "\n")
    cat("Model: independent\n")
    cat("Number of clustered individuals: ", ownNum2Str(ninds), "\n")
    cat("Number of groups in optimal partition: ", ownNum2Str(npops), "\n")
    cat("Log(marginal likelihood) of optimal partition: ", ownNum2Str(logml), "\n")
    cat(" ")
  }
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
  if (verbose) cat("Best Partition: ")
  if (fid != -1) {
    append(fid, c("Best Partition: ", "\n"))
  }
  for (m in 1:cluster_count) {
    indsInM <- matlab2r::find(PARTITION == m)
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
        text <- c(text, ",", as.character(indsInM[k]))
      }
    }
    text <- c(text, "}\n")
    while (length(text) > 58) {
      # Take one line and display it.
      new_line <- takeLine(text, 58)
      text <- (length(new_line) + 1):length(text)
      if (verbose) cat(new_line)
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
    if (any(text != "")) {
      if (verbose) cat(text)
      if (fid != -1) {
        append(fid, text)
        append(fid, "\n")
      }
    }
  }

  if (npops > 1) {
    if (verbose) {
      cat("\n")
      cat("\n")
      cat(
        "Changes in log(marginal likelihood)",
        " if indvidual i is moved to group j:\n"
      )
    }
    if (fid != -1) {
      append(fid, " ")
      append(fid, "\n")
      append(fid, " ")
      append(fid, "\n")
      append(
        fid,
        c(
          "Changes in log(marginal likelihood)",
          "if indvidual i is moved to group j:\n"
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
      maxSize <- base::max(nameSizes)
      maxSize <- base::max(maxSize, 5)
      erotus <- maxSize - 5
      alku <- blanks(erotus)
      ekarivi <- c(alku, "  ind", blanks(6 + erotus))
    } else {
      ekarivi <- "  ind      "
    }
    for (i in 1:cluster_count) {
      ekarivi <- c(ekarivi, ownNum2Str(i), blanks(8 - floor(log10(i))))
    }
    if (verbose) cat(ekarivi)
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
        rivi <- c(blanks(maxSize - length(nimi)), nimi, ":\n")
      } else {
        rivi <- c("\n", blanks(4 - floor(log10(ind))), ownNum2Str(ind), ":\n")
      }
      for (j in 1:npops) {
        rivi <- c(rivi, "  ", logml2String(omaRound(muutokset[j])))
      }
      if (verbose) cat(rivi)
      if (fid != -1) {
        append(fid, rivi)
        append(fid, "\n")
      }
    }
    if (verbose) cat("\n\nKL-divergence matrix in PHYLIP format:\n")

    dist_mat <- zeros(npops, npops)
    if (fid != -1) {
      append(fid, " ")
      append(fid, " ")
      append(fid, "KL-divergence matrix in PHYLIP format:")
      append(fid, "\n")
    }

    COUNTS <- COUNTS[seq_len(nrow(adjprior)), seq_len(ncol(adjprior)), , drop = FALSE]
    maxnoalle <- size(COUNTS, 1)
    nloci <- size(COUNTS, 2)
    d <- zeros(maxnoalle, nloci, npops)
    prior <- adjprior
    prior[matlab2r::find(prior == 1)] <- 0

    # Loci in which only one allele was detected.
    nollia <- matlab2r::find(all(prior == 0))

    prior[1, nollia] <- 1
    for (pop1 in 1:npops) {
      squeezed_COUNTS_prior <- squeeze(COUNTS[, , pop1]) + prior
      d[, , pop1] <- squeezed_COUNTS_prior / sum(squeezed_COUNTS_prior)
    }
    ekarivi <- as.character(npops)
    if (verbose) cat(ekarivi)
    if (fid != -1) {
      append(fid, ekarivi)
      append(fid, "\n")
    }

    for (pop1 in 1:npops) {
      for (pop2 in seq_len(pop1 - 1)) {
        dist1 <- d[, , pop1]
        dist2 <- d[, , pop2]
        div12 <- sum(
          sum(dist1 * base::log2((dist1 + 10^-10) / (dist2 + 10^-10)))
        ) / nloci
        div21 <- sum(
          sum(dist2 * base::log2((dist2 + 10^-10) / (dist1 + 10^-10)))
        ) / nloci
        div <- (div12 + div21) / 2
        dist_mat[pop1, pop2] <- div
      }
    }


    dist_mat <- dist_mat + t(dist_mat) # make it symmetric
    for (pop1 in 1:npops) {
      rivi <- c("\nCluster_", as.character(pop1), "\n")
      for (pop2 in 1:npops) {
        rivi <- c(rivi, kldiv2str(dist_mat[pop1, pop2]))
      }
      if (verbose) cat(rivi)
      if (fid != -1) {
        append(fid, rivi)
        append(fid, "\n")
      }
    }
  }
  if (verbose) {
    cat(
      "\n\nList of sizes of 10 best visited partitions",
      "and corresponding log(ml) values\n"
    )
  }

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
  partitionSummary <- partitionSummary[matlab2r::find(partitionSummary[, 2] > -1e49), ]
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
    if (verbose) cat(line)
    if (fid != -1) {
      append(fid, line)
      append(fid, "\n")
    }
  }

  if (!fixedK) {
    if (verbose) cat("\n\nProbabilities for number of clusters\n")
    if (fid != -1) {
      append(fid, " ")
      append(fid, "\n")
      append(fid, " ")
      append(fid, "\n")
      append(fid, "Probabilities for number of clusters")
      append(fid, "\n")
    }

    npopsTaulu <- unique(partitionSummary[, 1])
    len <- length(npopsTaulu)
    probs <- zeros(len, 1)
    partitionSummary[, 2] <- partitionSummary[, 2] -
      base::max(partitionSummary[, 2])
    sumtn <- sum(exp(partitionSummary[, 2]))
    for (i in 1:len) {
      npopstn <- sum(
        exp(
          partitionSummary[matlab2r::find(
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
        if (verbose) cat(line, "\n")
        if (fid != -1) {
          append(fid, line)
          append(fid, "\n")
        }
      }
    }
  }
  # Closing sink(s)
  while (sink.number() > 0L) {
    sink()
  }
  return(changesInLogml)
}
