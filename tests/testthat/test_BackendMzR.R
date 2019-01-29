test_that(".spectra_from_file_mzR works", {
    spd <- fData(sciex)
    expect_error(.spectra_from_file_mzR())
    res <- .spectra_from_file_mzR(fileNames(sciex)[1],
                                  spd[spd$fileIdx == 1, ])
    expect_equal(res, spectra(filterFile(sciex, 1)))
    expect_true(length(res) == 931)
    res <- .spectra_from_file_mzR(fileNames(sciex)[1],
                                  spd[13, ])
    expect_equal(res[[1]], sciex[[13]])
    expect_true(length(mz(res[[1]])) == 1650)
})

test_that("backendSpectrapply,BackendMzR works", {
    sciex_spectra <- spectra(sciex) # TODO: Change this to sciex_inmem
    spd <- sciex_mzr@spectraData
    be <- BackendMzR()
    be <- backendInitialize(be, files = sf)
    expect_true(validObject(be))
    ## Get all spectra from a single file
    res <- backendSpectrapply(be, spectraData = spd[spd$fileIdx == 2, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 2, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 2])

    ## Get some spectra from both file s
    idx <- c(4, 123, 432, 1110, 1234)
    res <- backendSpectrapply(be, spectraData = spd[idx, ])
    expect_equal(names(res), rownames(spd)[idx])
    expect_equal(res, sciex_spectra[idx])

    ## With FUN.
    res <- backendSpectrapply(be, spectraData = spd[idx, ],
                              FUN = function(z) mean(mz(z)))
    expect_equal(names(res), rownames(spd)[idx])
    expect_true(is.list(res))
    expect_true(all(vapply(res, is.numeric, logical(1))))
    expect_equal(res, lapply(sciex_spectra[idx], function(z) mean(mz(z))))

    ## With processingStep
    be@processingQueue <- list(ProcessingStep("mz"))
    res <- backendSpectrapply(be, spectraData = spd[idx, ])
    expect_equal(names(res), rownames(spd)[idx])
    expect_true(is.list(res))
    expect_true(all(vapply(res, is.numeric, logical(1))))
    expect_equal(res, lapply(sciex_spectra[idx], mz))

    ## With processingStep and FUN
    res <- backendSpectrapply(be, spectraData = spd[idx, ],
                              FUN = "mean")
    expect_equal(names(res), rownames(spd)[idx])
    expect_true(is.list(res))
    expect_true(all(vapply(res, is.numeric, logical(1))))
    expect_equal(res, lapply(sciex_spectra[idx], function(z) mean(mz(z))))

    ## Errors
    be@files <- c("a", "b")
    expect_error(backendSpectrapply(be, spectraData = spd[idx, ]))
    be@files <- "d"
    expect_error(backendSpectrapply(be, spectraData = spd[idx, ]))
})

test_that("backendReadSpectra,BackendMzR works", {
    sciex_spectra <- spectra(sciex) # TODO: Change this to sciex_inmem
    spd <- sciex_inmem@spectraData
    be <- BackendMzR()
    be <- backendInitialize(be, sf)

    ## Get spectra from a single file
    res <- backendReadSpectra(be, spd[spd$fileIdx == 1, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 1, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 1])
    res <- backendReadSpectra(be, spd[spd$fileIdx == 2, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 2, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 2])
})

test_that("backendAddProcessing,Backend default works", {
    be <- sciex_mzr@backend
    expect_true(length(be@processingQueue) == 0)
    be <- MSnbase:::backendAddProcessing(be, procStep = ProcessingStep("mz"))
    expect_true(length(be@processingQueue) == 1)
    tmp <- sciex_mzr
    tmp@backend <- be
    res <- spectrapply(tmp)
    expect_true(all(vapply(res, is.numeric, logical(1))))
})
