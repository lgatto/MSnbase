context("MSnProcess class")

test_that("MSnbase version", {
  p <- new("MSnProcess")
  expect_true(validObject(p))
  expect_true(p@MSnbaseVersion==as.character(packageVersion("MSnbase")))
})
