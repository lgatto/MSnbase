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

test_that("backendApplyProcessingQueue,BackendMzR works", {
    be <- sciex_mzr@backend
    spd <- sciex_mzr@spectraData[c(13, 15, 33, 113, 117, 167), ]
    the_q <- list(ProcessingStep("removePeaks", list(t = 5000)))
    v <- isMSnbaseVerbose()
    setMSnbaseVerbose(TRUE)
    expect_message(res <- backendApplyProcessingQueue(be, spd, queue = the_q))
    expect_equal(res@modCount, c(0L, 0L))
    setMSnbaseVerbose(v)
})
