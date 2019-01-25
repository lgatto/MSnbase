test_that("MSnExperiment validator works", {
    library(MSnbase)
    library(testthat)
    tst <- new("MSnExperiment")
    expect_error(validObject(tst))
    tst <- new("MSnExperiment", backend = BackendMzR())
    expect_true(validObject(tst))
    tst@files <- "a"
    expect_true(is.character(.validMSnExperiment(tst)))
    expect_error(validObject(tst))
    tst@sampleData <- DataFrame(sampleIdx = 1)
    expect_true(.validMSnExperiment(tst))
})

test_that("readMSnExperiment works", {
    sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    res <- readMSnExperiment(sf)
    expect_true(validObject(res))
    expect_true(length(unique(res@spectraData$fileIdx)) == 2)
    expect_error(readMSnExperiment(sf, sampleData = DataFrame(a = 3)))
})
