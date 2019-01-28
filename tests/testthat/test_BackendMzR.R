test_that(".spectra_from_file_mzR works", {
    spd <- fData(sciex)
    expect_error(.spectra_from_file_mzR())
    res <- MSnbase:::.spectra_from_file_mzR(fileNames(sciex)[1],
                                            spd[spd$fileIdx == 1, ])
    expect_equal(res, spectra(filterFile(sciex, 1)))
    expect_true(length(res) == 931)
    res <- MSnbase:::.spectra_from_file_mzR(fileNames(sciex)[1],
                                            spd[13, ])
    expect_equal(res[[1]], sciex[[13]])
    expect_true(length(mz(res[[1]])) == 1650)
})

test_that("backendReadSpectra works", {
    sciex_spectra <- spectra(sciex) # Use the in-mem backend later
    sciex_me <- readMSnExperiment(sf)
    spd <- sciex_me@spectraData
    fls <- sciex_me@files
    be <- BackendMzR()
    ## Errors
    expect_error(MSnbase:::backendReadSpectra(be, file = c("a"), spd))
    expect_error(MSnbase:::backendReadSpectra(be, file = c("a"),
                                              spd[spd$fileIdx == 1, ]))

    ## Get spectra from a single file
    res <- MSnbase:::backendReadSpectra(be, file = fls[1],
                                        spd[spd$fileIdx == 1, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 1, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 1])
    res <- MSnbase:::backendReadSpectra(be, file = fls[2],
                                        spd[spd$fileIdx == 2, ])
    expect_equal(names(res), rownames(spd[spd$fileIdx == 2, ]))
    expect_true(all(vapply(res, is, "Spectrum", FUN.VALUE = logical(1))))
    expect_equal(res, sciex_spectra[spd$fileIdx == 2])

    ## Get spectra from multiple files
    res <- MSnbase:::backendReadSpectra(be, file = fls, spd)
    expect_equal(res, sciex_spectra)
})
