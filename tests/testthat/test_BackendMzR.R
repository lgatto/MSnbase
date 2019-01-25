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
