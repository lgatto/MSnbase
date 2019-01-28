test_that("MSnExperiment validator works", {
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

test_that("spectrapply,MSnExperiment works", {
    sciex_me <- readMSnExperiment(sf)
    sciex_spctra <- spectra(sciex)

    res <- spectrapply(sciex_me)
    expect_equal(names(res), rownames(sciex_me@spectraData))
    expect_true(all(vapply(res, is, logical(1), "Spectrum")))
    expect_true(all(vapply(res, validObject, logical(1))))
    expect_equal(sciex_spctra, res)

    res_2 <- spectrapply(sciex_me, FUN = function(z) mean(mz(z)))
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))

    sciex_me@processingQueue <- list(ProcessingStep("mz"))
    res_2 <- spectrapply(sciex_me, FUN = "mean")
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))

    sciex_me@processingQueue <- list()
    res <- spectrapply(sciex_inmem)
    expect_equal(names(res), rownames(sciex_me@spectraData))
    expect_true(all(vapply(res, is, logical(1), "Spectrum")))
    expect_true(all(vapply(res, validObject, logical(1))))
    expect_equal(sciex_spctra, res)
})
