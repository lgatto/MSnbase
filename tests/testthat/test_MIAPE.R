context("MIAPE class")

data(itraqdata)

test_that("MIAPE validity", {
    expect_true(validObject(new("MIAPE")))
    expect_null(show(experimentData(itraqdata)))
    itraqdata@experimentData@instrumentModel <- "test"
    itraqdata@experimentData@other <- list(Note1 = "First note")
    expect_null(show(experimentData(itraqdata)))
    expect_null(msInfo(itraqdata))
})

test_that("misc metadata", {
    expect_identical(analyzerDetails(itraqdata),
                     analyserDetails(itraqdata))
    expect_identical(exptitle(itraqdata), "Example 'MSnExp' data set")
    expect_identical(expemail(itraqdata), "lg390@cam.ac.uk")
})
