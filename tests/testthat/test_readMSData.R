test_that("Empty data", {
    f <- msdata::proteomics(full.names = TRUE,
                            pattern = "MRM-standmix-5.mzML")
    expect_warning(x <- readMSData(f, mode = "onDisk"))
    expect_identical(length(x), 0L)
    expect_true(inherits(x, "OnDiskMSnExp"))
    expect_warning(x <- readMSData(f, mode = "inMemory"))
    expect_identical(length(x), 0L)
    expect_true(inherits(x, "MSnExp"))
})

test_that("One empty data file", {
    f1 <- msdata::proteomics(full.names = TRUE,
                             pattern = "MRM-standmix-5.mzML")
    f2 <- msdata::proteomics(full.names = TRUE,
                             pattern = "MS3TMT11.mzML")
    expect_warning(x <- readMSData(c(f1, f2), mode = "onDisk"))
    expect_identical(fileNames(x), f2)
})
