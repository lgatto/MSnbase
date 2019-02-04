sciex_h5 <- readMSnExperiment(sf, backend = BackendHdf5(),
                              path = paste0(tempdir(), "/sciex_h5/"))

test_that("BackendHdf5 validators work", {
    expect_true(is.null(.valid.BackendHdf5.checksums(5, 5)))
    expect_true(is.character(.valid.BackendHdf5.checksums(3, 1:2)))
    expect_true(is.null(.valid.BackendHdf5.h5files(5, 5)))
    expect_true(is.character(.valid.BackendHdf5.h5files(3, 1:2)))
    expect_true(is.character(.valid.BackendHdf5.h5files(c(1, 1), c(1, 2))))
    tmpf <- tempfile()
    write(4, tmpf)
    expect_true(is.null(.valid.BackendHdf5.h5files.exist(tmpf)))
    expect_true(is.character(.valid.BackendHdf5.h5files.exist("4")))
})

test_that("BackendHdf5 works", {
    tst <- BackendHdf5()
    expect_true(is(tst, "BackendHdf5"))
    expect_true(validObject(tst))
    tst@checksums <- "r"
    expect_error(validObject(tst))
})

test_that("backendInitialize and backendImportData,BackendHdf5 work", {
    tst <- BackendHdf5()
    dr <- tempdir()
    res <- MSnbase:::backendInitialize(tst, files = sf, path = dr)
    expect_equal(length(res@h5files), 2)
    expect_equal(length(res@checksums), 2)
    expect_true(all(file.exists(res@h5files)))
    expect_error(MSnbase:::backendInitialize(tst, files = sf, path = dr))
    expect_true(validObject(res))
    expect_true(all(res@checksums == c("", "")))
    show(res)

    ## backendImportData
    spdata <- sciex_mzr@spectraData
    res <- MSnbase:::backendImportData(res, spdata)
    cont <- rhdf5::h5ls(res@h5files[1])
    expect_equal(sum(spdata$fileIdx == 1) + 3, nrow(cont))
    cont <- rhdf5::h5ls(res@h5files[2])
    expect_equal(sum(spdata$fileIdx == 2) + 3, nrow(cont))
})

test_that(".serialize_msfile_to_hdf5 works", {
    h5file <- tempfile()
    h5 <- rhdf5::H5Fcreate(h5file)
    rhdf5::h5createGroup(h5, "spectra")
    rhdf5::h5createGroup(h5, "md5")
    rhdf5::H5Fclose(h5)
    md5 <- MSnbase:::.serialize_msfile_to_hdf5(fileNames(sciex)[1], h5file)
    cont <- rhdf5::h5ls(h5file)
    expect_equal(nrow(cont), sum(fromFile(sciex) == 1) + 3)
})

test_that(".h5_read_bare works", {
    expect_error(MSnbase:::.h5_read_bare())
    expect_error(MSnbase:::.h5_read_bare("5"))

    fid <- .Call("_H5Fopen", sciex_h5@backend@h5files[1], 0L, PACKAGE = "rhdf5")
    res <- MSnbase:::.h5_read_bare(fid, "/spectra/3")
    .Call("_H5Fclose", fid, PACKAGE = "rhdf5")
    res_2 <- rhdf5::h5read(sciex_h5@backend@h5files[1],
                           name = "/spectra/3")
    expect_equal(res, res_2)
    expect_equal(res[, 1], mz(sciex[[3]]))
    expect_equal(res[, 2], intensity(sciex[[3]]))
})

test_that(".h5_read_spectra works", {
    ## fd <- fData(h5_sciex)
    ## res_h5 <- MSnbase:::.h5_read_spectra(fd[fd$fileIdx == 1, , drop = FALSE],
    ##                                      h5_sciex@h5files[1],
    ##                                      h5_sciex@checksums[1],
    ##                                      fileNames(h5_sciex)[1])
    ## res_od <- spectra(sciex)[fromFile(sciex) == 1]
    ## expect_equal(res_h5, res_od)
    ## idx <- c(34, 65, 234, 453, 488)
    ## res_h5 <- MSnbase:::.h5_read_spectra(fd[idx, ], h5_sciex@h5files[1],
    ##                                      h5_sciex@checksums[1],
    ##                                      fileNames(h5_sciex)[1])
    ## expect_equal(res_h5, res_od[idx])
    ## ## MS1 & 2 data
    ## res_h5 <- MSnbase:::.h5_read_spectra(fData(h5_tmt), h5_tmt@h5files,
    ##                                      h5_tmt@checksums,
    ##                                      fileNames(h5_tmt))
    ## res_od <- spectra(tmt_erwinia_on_disk)
    ## expect_equal(res_h5, res_od)
    ## expect_equal(res_od[123], MSnbase:::.h5_read_spectra(fData(h5_tmt)[123, ],
    ##                                                      h5_tmt@h5files,
    ##                                                      h5_tmt@checksums,
    ##                                                      fileNames(h5_tmt)))
})
