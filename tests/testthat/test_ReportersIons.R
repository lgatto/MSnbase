context("ReporterIons class")

test_that("ReporterIons validity", {
  expect_true(validObject(iTRAQ4))
  expect_true(validObject(iTRAQ5))
  expect_true(validObject(iTRAQ8))
  expect_true(validObject(iTRAQ9))
  expect_true(validObject(TMT6))
  expect_true(validObject(TMT7))
  expect_true(validObject(TMT10))
  expect_null(show(TMT10))
})

test_that("Reporter data", {
    expect_identical(reporterNames(iTRAQ4),
                     paste("iTRAQ4", 114:117, sep = "."))
    expect_identical(reporterColours(iTRAQ4),
                     c("red", "green", "blue", "yellow"))
    expect_identical(reporterColours(iTRAQ4),
                     reporterColors(iTRAQ4))
    reporterNames(iTRAQ4) <- as.character(1:4)
    expect_identical(reporterNames(iTRAQ4), as.character(1:4))
})
