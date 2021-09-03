context("Auxiliary functions to greedyMix")

# Defining the relative path to current inst -----------------------------------
if (interactive()) {
	path_inst <- "../../inst/ext"
} else {
	path_inst <- system.file("ext", "", package="rBAPS")
}

# Reading datasets -------------------------------------------------------------
baps_diploid <- read.delim(
	file = paste(path_inst, "BAPS_format_clustering_diploid.txt", sep="/"),
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

context("Opening files on greedyMix")

df_fasta <- greedyMix(
	data   = file.path(path_inst, "FASTA_clustering_haploid.fasta"),
	format = "FASTA"
)
df_vcf <- greedyMix(
	data    = file.path(path_inst, "vcf_example.vcf"),
	format  = "VCF",
	verbose = FALSE
)
df_bam <- greedyMix(
	data    = file.path(path_inst, "bam_example.bam"),
	format  = "BAM",
)
# TODO #19: add example reading Genpop
test_that("Files are imported correctly", {
	expect_equal(dim(df_fasta), c(5, 99))
	expect_equal(dim(df_vcf), c(variants = 2, fix_cols = 8, gt_cols = 3))
	expect_error(
		greedyMix(
			data    = paste(path_inst, "sam_example.sam", sep="/"),
			format  = "SAM",
		)
	)
	expect_equal(length(df_bam[[1]]), 13)
})

context("Linkage")

test_that("Linkages are properly calculated", {
	Y <- c(0.5, 0.3, 0.6, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3, 0.5)
	expect_equal(
		object   = linkage(Y),
		expected = matrix(c(2, 1, 7, 8, 4, 3, 5, 6, .2, .3, .3, .6), ncol=3)
	)
})
