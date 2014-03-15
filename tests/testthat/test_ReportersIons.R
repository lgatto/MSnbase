context("ReporterIons class")

test_that("ReporterIons validity", {
  expect_true(validObject(iTRAQ4))
  expect_true(validObject(iTRAQ5))
  expect_true(validObject(iTRAQ8))
  expect_true(validObject(iTRAQ9))
  expect_true(validObject(TMT6))
  expect_true(validObject(TMT7))
})
