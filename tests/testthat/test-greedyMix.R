context("Auxiliary functions to greedyMix")

baps_diploid <- read.delim(
	"inst/ext/ExamplesDataFormatting/Example data in BAPS format for clustering of diploid individuals.txt",
	sep = " ",
	header = FALSE
)

handleData(baps_diploid)$newData

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

# TODO: needs #12 to be fixed before this can be done without user intervention
greedyMix(
	tietue = "inst/ext/ExamplesDataFormatting/Example data in BAPS format for clustering of diploid individuals.txt",
	format = "BAPS",
	savePreProcessed = FALSE
) # Upper bounds 100 100

context("Linkage")

test_that("Linkages are properly calculated", {
	Y <- c(0.5, 0.3, 0.6, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3, 0.5)
	expect_equal(
		object   = linkage(Y),
		expected = matrix(c(2, 1, 7, 8, 4, 3, 5, 6, .2, .3, .3, .6), ncol=3)
	)
})