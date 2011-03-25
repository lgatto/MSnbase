context("MSnSet class")

test_that("MSnSet validity", {
  expect_true(validObject(new("MSnSet")))
})
