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
    be <- MSnbase:::backendInitialize(be, sf)

    ## Get spectra from a single file
    res <- MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 1, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 1, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 1])
    res <- MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 2, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 2, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 2])
})
