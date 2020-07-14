context("Basic Matlab functions")

test_that("rand works properly", {
	expect_equal(dim(rand()), c(1, 1))
	expect_equal(dim(rand(1, 2)), c(1, 2))
	expect_equal(dim(rand(3, 2)), c(3, 2))
})

test_that("repmat works properly", {
	mx0 <- c(1:4) # when converted to matrix, results in a column vector
	mx1 <- matrix(5:8)
	mx2 <- matrix(0:-3, 2)
	expect_error(repmat(mx0))
	expect_equal(repmat(mx0, 1), t(as.matrix(mx0)))
	expect_equal(
		object = repmat(mx0, 2),
		expected = unname(cbind(rbind(mx0, mx0), rbind(mx0, mx0)))
	)
	expect_equal(
		object = repmat(mx1, 2),
		expected = unname(cbind(rbind(mx1, mx1), rbind(mx1, mx1)))
	)
	expect_equal(
		object = repmat(mx2, c(2, 3)),
		expected = cbind(rbind(mx2, mx2), rbind(mx2, mx2), rbind(mx2, mx2))
	)
	expect_equal(
		object = repmat(mx2, c(4, 1)),
		expected = rbind(mx2, mx2, mx2, mx2)
	)
	expect_equal(
		object = repmat(mx2, c(1, 1, 2)),
		expected = array(mx2, c(2, 2, 2))
	)
	expect_equal(repmat(1:2, 3), matrix(rep(1:2, 9), 3, 6, byrow=TRUE))
	expect_equal(repmat(10, c(3, 2)), matrix(10, 3, 2))
})

test_that("zeros and ones work as expected", {
	expect_equal(zeros(1), matrix(0, 1))
	expect_equal(zeros(2), matrix(0, 2, 2))
	expect_equal(zeros(2, 1), matrix(0, 2, 1))
	expect_equal(zeros(1, 10), matrix(0, 1, 10))
	expect_equal(zeros(3, 2, 4), array(0, c(3, 2, 4)))
	expect_equal(ones(8), matrix(1, 8, 8))
	expect_equal(ones(5, 2), matrix(1, 5, 2))
	expect_equal(ones(2, 100), matrix(1, 2, 100))
	expect_equal(ones(3, 2, 4, 2), array(1, c(3, 2, 4, 2)))
})

test_that("times works as expected", {
	expect_equal(times(9, 6), as.matrix(54))
	expect_equal(times(9, c(.8, 9)), as.matrix(c(7.2, 81)))
	expect_equal(times(c(.8, 9), 5), as.matrix(c(4, 45)))
	expect_equal(times(matrix(1:4, 2), 5), matrix(c(5, 10, 15, 20), 2))
	expect_equal(times(5, matrix(1:4, 2)), matrix(c(5, 10, 15, 20), 2))
	expect_equal(times(matrix(1:4, 2), c(10, 3)), matrix(c(10, 6, 30, 12), 2))
	expect_equal(
		object = times(matrix(1:4, 2), matrix(c(10, 3), 1)),
		expected = matrix(c(10, 20, 9, 12), 2)
	)
	expect_equal(times(c(10, 3), matrix(1:4, 2)), matrix(c(10, 6, 30, 12), 2))
	expect_equal(
		object = times(matrix(c(10, -5, 3, 9), 2), matrix(1:4, 2)),
		expected = matrix(c(10, -10, 9, 36), 2)
	)
	expect_equal(
		object = times(matrix(c(-1.6, 5), 1), c(8, 1)),
		expected = matrix(c(-12.8, -1.6, 40, 5), 2)
	)
})

test_that("colon works as expected (hee hee)", {
	expect_equal(colon(1, 4), 1:4)
	expect_length(colon(4, 1), 0)
})

test_that("size works as on MATLAB", {
	sk <- 10
	vk <- 1:4
	mx <- matrix(1:6, 2)
	ra <- array(1:24, c(2, 3, 4))
	expect_equal(size(sk), 1)
	expect_equal(size(vk), c(1, 4))
	expect_equal(size(mx), c(2, 3))
	expect_equal(size(ra), c(2, 3, 4))
	expect_equal(size(sk, 199), 1)
	expect_equal(size(vk, 199), 1)
	expect_equal(size(mx, 199), 1)
	expect_equal(size(ra, 199), 1)
	expect_equal(size(vk, 2), 4)
	expect_equal(size(mx, 2), 3)
	expect_equal(size(ra, 2), 3)
	expect_equal(size(ra, 3), 4)
})

test_that("reshape reshapes properly", {
	mx <- matrix(1:4, 2)
	ra <- array(1:12, c(2, 3, 2))
	expect_equal(reshape(mx, c(1, 4)), matrix(1:4, 1))
	expect_equal(reshape(mx, c(2, 2)), mx)
	expect_equal(reshape(mx, c(1, 1, 4)), array(mx, c(1, 1, 4)))
	expect_error(reshape(mx, c(1, 2, 3)))
	expect_error(reshape(ra, c(1, 2, 3)))
	expect_equal(reshape(ra, c(3, 2, 2)), array(ra, c(3, 2, 2)))
})

test_that("isfield works as on Matlab", {
	S <- list()
	S$x <- rnorm(100)
	S$y <- sin(S$x)
	S$title <- "y = sin(x)"
	expect_true(isfield(S, "title"))
	expect_equivalent(
		object = isfield(S, c("x", "y", "z", "title", "error")),
		expected = c(TRUE, TRUE, FALSE, TRUE, FALSE)
	)
})

test_that("strcmp works as expected", {
	yes <- 'Yes'
	no <- 'No'
	ja <- 'Yes'
	expect_false(strcmp(yes, no))
	expect_true(strcmp(yes, ja))
	s1 <- 'upon'
	s2 <- matrix(c('Once', 'upon', 'a', 'time'), 2, byrow=TRUE)
	s3 <- c('Once', 'upon', 'a', 'time')
	s4 <- matrix(c("A", "bc", "def", "G"), 2, byrow=TRUE)
	s5 <- matrix(c("B", "c", "def", "G"), 2, byrow=TRUE)
	expect_equal(strcmp(s1, s2), matrix(c(FALSE, FALSE, TRUE, FALSE), 2))
	expect_equivalent(strcmp(s1, s3), c(FALSE, TRUE, FALSE, FALSE))
	expect_error(strcmp(s2, s3))
	expect_equal(strcmp(s4, s5), matrix(c(FALSE, TRUE, FALSE, TRUE), 2))
})

test_that("isempty works as expected", {
	A <- array(dim=c(0, 2, 2))
	B <- matrix(rep(NA, 4), 2)
	C <- matrix(rep(0, 4), 2)
	cat1 <- as.factor(c(NA, NA))
	cat2 <- as.factor(c())
	str1 <- matrix(rep("", 3))
	expect_true(isempty(A))
	expect_false(isempty(B))
	expect_false(isempty(C))
	expect_false(isempty(cat1))
	expect_true(isempty(cat2))
	expect_false(isempty(str1))
})

test_that("find works as expected", {
	X <- matrix(c(1, 0, 2, 0, 1, 1, 0, 0, 4), 3, byrow=TRUE)
	Y <- seq(1, 19, 2)
	expect_equal(find(X), c(1, 5, 7, 8, 9))
	expect_equal(find(!X), c(2, 3, 4, 6))
	expect_equal(find(Y == 13), 7)
})

test_that("sortrows works as expected", {
	 mx <- matrix(c(3, 2, 2, 1, 1, 10, 0, pi), 4)
	 expect_equal(sortrows(mx), matrix(c(1, 2, 2, 3, pi, 10, 0, 1), 4))
	 expect_equal(sortrows(mx, 2), matrix(c(2, 3, 1, 2, 0, 1, pi, 10), 4))
	 expect_equal(sortrows(mx, 1:2), mx[order(mx[, 1], mx[, 2]), ])
})

test_that("cell works as expected", {
	expect_equal(cell(0), array(dim = c(0, 0)))
	expect_equal(cell(1), array(dim = c(1, 1)))
	expect_equal(cell(2), array(dim = c(2, 2)))
	expect_equal(cell(3, 4), array(dim = c(3, 4)))
	expect_equal(cell(5, 7, 6), array(dim = c(5, 7, 6)))
})

test_that("blanks works as expected", {
	expect_warning(blanks(-1))
	expect_equal(suppressWarnings(blanks(-1)), "")
	expect_equal(blanks(0), "")
	expect_equal(blanks(1), " ")
	expect_equal(blanks(10), "          ")
})