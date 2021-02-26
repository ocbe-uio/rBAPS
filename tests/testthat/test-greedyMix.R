context("Opening files on greedyMix")

# TODO: needs #12 to be fixed before this can be done without user intervention
greedyMix(
	tietue = "inst/ext/ExamplesDataFormatting/Example data in BAPS format for clustering of diploid individuals.txt",
	format = "GenePop",
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