context("MIAPE class")

test_that("MIAPE validity", {
    expect_true(validObject(new("MIAPE")))
    data(itraqdata)
    expect_null(show(experimentData(itraqdata)))
    expect_null(msInfo(itraqdata))
})
