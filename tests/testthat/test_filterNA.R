context("filterNA method")

test_that("filterNA for matrix", {
  M <- matrix(rnorm(20, 1), 5, 4)
  M[2, 1] <- NA
  M[3, 1:2] <- NA
  M[4, 1:3] <- NA
  M[5, ] <- NA
  rownames(M) <- LETTERS[1:5]
  expect_equal(rownames(filterNA(M, 0/4)), LETTERS[1])
  expect_equal(rownames(filterNA(M, 1/4)), LETTERS[1:2])
  expect_equal(rownames(filterNA(M, 2/4)), LETTERS[1:3])
  expect_equal(rownames(filterNA(M, 3/4)), LETTERS[1:4])
  expect_equal(rownames(filterNA(M, 4/4)), LETTERS[1:5])
  M[1, 1] <- NA
  expect_equal(dim(filterNA(M, 0/4)), c(0, 4))
})

test_that("filterNA for MSnSet", {
  xx <- quantify(itraqdata, reporters = iTRAQ4, method = "max", verbose=FALSE)
  exprs(xx)[sample(prod(dim(xx)), 120)] <- NA
  xx2 <- filterNA(xx, 1/4)
  expect_true(nrow(xx) > nrow(xx2))
})

test_that("filterNA with patter", {
  
  M <- matrix(rnorm(40), ncol = 4)
  rownames(M) <- 1:10
  M[1:5, 1] <- NA
  M[c(1, 3, 5), 4] <- NA
  M[c(1, 5, 10), 2] <- NA
  res <- filterNA(M, pattern = "0110")
  expect_equal(rownames(res), as.character(c(2, 3, 4, 6, 7, 8, 9)))
  res <- filterNA(M, pattern = "0111")
  expect_equal(rownames(res), as.character(c(2, 4, 6, 7, 8, 9)))
  res0 <- filterNA(M, pattern = "1111")
  res1 <- filterNA(M, pNA = 0)
  expect_identical(res0, res1)
})
