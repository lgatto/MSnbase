context("Updating objects")

test_that("update MSnExp version 0.3.0 -> 0.3.1", {
    ## msx0 was created before Spectrum2 object got a polarity slot
    f0 <- dir(system.file(package = "MSnbase", dir = "extdata"),
              full.name = TRUE, pattern = "msx0.rda")
    load(f0) ## msx0
    expect_error(validObject(msx0))
    ## See issue #166 why this is disabled
    ## expect_true(validObject(updateObject(msx0)))
})
