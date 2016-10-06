test_that("centroided accessor with/without na.fail", {
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    x <- readMSData(f)
    x2 <- readMSData2(f)
    ## in-mem, single spectrum
    expect_true(is.na(centroided(x[[1]])))
    expect_error(centroided(x[[1]], na.fail = TRUE))
    ## in-mem, experiment
    expect_true(all(is.na(centroided(x))))
    expect_error(centroided(x, na.fail = TRUE))
    ## on-disk, single spectrum
    expect_true(is.na(centroided(x2[[1]])))
    expect_error(centroided(x2[[1]], na.fail = TRUE))
    ## on-disk, experiment
    expect_true(all(is.na(centroided(x2))))
    expect_error(centroided(x2, na.fail = TRUE))
    fData(x2)$centroided[1:(length(x2) - 1)] <- TRUE
    expect_error(centroided(x2, na.fail = TRUE))
})
