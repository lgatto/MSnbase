test_that("MSnExperiment validator works", {
    tst <- new("MSnExperiment")
    expect_error(validObject(tst))
    tst <- new("MSnExperiment", backend = BackendMzR())
    expect_true(validObject(tst))
    tst@sampleData <- DataFrame(sampleIdx = 1)
    res <- MSnbase:::.validMSnExperiment(tst)
    expect_true(is.character(res))
    expect_error(validObject(tst))
})

test_that("readMSnExperiment works", {
    sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    res <- readMSnExperiment(sf)
    expect_true(validObject(res))
    expect_true(length(unique(res@spectraData$fileIdx)) == 2)
    expect_error(readMSnExperiment(sf, sampleData = DataFrame(a = 3)))
})

test_that("spectrapply,MSnExperiment works", {
    sciex_spctra <- spectra(sciex)

    ## BackendMzR
    res <- spectrapply(sciex_mzr)
    expect_equal(names(res), rownames(sciex_mzr@spectraData))
    expect_true(all(vapply(res, is, logical(1), "Spectrum")))
    expect_true(all(vapply(res, validObject, logical(1))))
    expect_equal(sciex_spctra, res)

    res_2 <- spectrapply(sciex_mzr, FUN = function(z) mean(mz(z)))
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))

    sciex_mzr@backend@processingQueue <- list(ProcessingStep("mz"))
    res_2 <- spectrapply(sciex_mzr, FUN = "mean")
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))
    sciex_mzr@backend@processingQueue <- list()

    ## Compare BackendMzR and BackendMemory
    expect_equal(spectrapply(sciex_inmem), spectrapply(sciex_mzr))
    expect_equal(spectrapply(sciex_inmem, FUN = mz),
                 spectrapply(sciex_mzr, FUN = mz))
})
