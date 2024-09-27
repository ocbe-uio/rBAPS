context("Auxiliary functions to greedyMix")

# Defining the relative path to current inst -----------------------------------
path_inst <- system.file("extdata", "", package = "rBAPS")

# Reading datasets -------------------------------------------------------------
baps_diploid <- read.delim(
  file = file.path(path_inst, "BAPS_clustering_diploid.txt"),
  sep = " ",
  header = FALSE
)

test_that("handleData works as expected", {
  data_obs <- handleData(baps_diploid)$newData
  data_exp <- matrix(
    c(
      -9, 1, 2, 1, 1, 1, 2, 1, 2, 2, 1,
      -9, 1, 1, 2, 2, 2, 1, 1, 1, 2, 1,
      3, 2, 2, 3, 2, -9, 3, 1, 2, 1, 2,
      2, 1, 2, 1, 2, -9, 1, 1, 1, 1, 2,
      3, 1, 1, 1, 2, 1, 1, 2, -9, 1, 3,
      3, 1, 2, 1, 1, 1, 2, 1, -9, 2, 3,
      1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 4,
      3, 2, 2, 3, 2, 2, 3, 1, 2, 1, 4,
      2, 1, 2, 1, -9, 1, 1, 1, 1, 1, 5,
      3, 1, 1, 1, -9, 1, 1, 2, 1, 1, 5
    ),
    nrow = 10, byrow = TRUE
  )
  colnames(data_exp) <- colnames(data_obs)
  expect_equal(data_obs, data_exp)
})

context("Processing files through greedyMix")

raw_fasta <- importFile(
  data   = file.path(path_inst, "FASTA_clustering_haploid.fasta"),
  format = "FASTA"
)
raw_vcf <- importFile(
  data    = file.path(path_inst, "vcf_example.vcf"),
  format  = "VCF",
  verbose = FALSE
)
raw_bam <- importFile(
  data    = file.path(path_inst, "bam_example.bam"),
  format  = "BAM",
)

test_that("Files are imported correctly", {
  expect_equal(dim(raw_fasta), c(5, 99))
  expect_equal(dim(raw_vcf), c(variants = 2, fix_cols = 8, gt_cols = 3))
  expect_error(
    importFile(
      data    = paste(path_inst, "sam_example.sam", sep = "/"),
      format  = "SAM",
    )
  )
  expect_equal(length(raw_bam[[1]]), 13)
})

test_that("greedyMix() fails successfully", {
  expect_error(greedyMix(file.path(path_inst, "vcf_example.vcf")))
  expect_error(greedyMix(file.path(path_inst, "bam_example.bam")))
})

test_that("greedyMix() works when it should", {
  tol <- 1e-4
  baps_file <- file.path(path_inst, "BAPS_clustering_diploid.txt")
  greedy_baps <- greedyMix(baps_file, "BAPS")
  expect_type(greedy_baps, "list")
  expect_length(greedy_baps, 10L)
  expect_equal(greedy_baps[['logml']], -78.10151, tolerance = tol)

  genepop_file <- file.path(path_inst, "GenePop.txt")
  greedy_genepop <- greedyMix(genepop_file, "GenePop")
  expect_type(greedy_genepop, "list")
  expect_length(greedy_genepop, 10L)
  expect_equal(greedy_genepop[['logml']], -10.76224, tolerance = tol)

  fasta_file <- file.path(path_inst, "FASTA_clustering_haploid.fasta")
})

context("Linkage")

test_that("Linkages are properly calculated", {
  Y <- c(0.5, 0.3, 0.6, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3, 0.5)
  expect_equal(
    object   = linkage(Y),
    expected = matrix(c(2, 1, 7, 8, 4, 3, 5, 6, .2, .3, .3, .6), ncol = 3)
  )
})
