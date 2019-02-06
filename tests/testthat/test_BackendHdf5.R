
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
    res <- backendInitialize(tst, files = sf, path = dr)
    expect_equal(length(res@h5files), 2)
    expect_equal(length(res@checksums), 2)
    expect_true(all(file.exists(res@h5files)))
    expect_error(backendInitialize(tst, files = sf, path = dr))
    expect_true(validObject(res))
    expect_true(all(res@checksums == c("", "")))
    show(res)

    ## backendImportData
    spdata <- sciex_mzr@spectraData
    res <- backendImportData(res, spdata)
    cont <- rhdf5::h5ls(res@h5files[1])
    expect_equal(sum(spdata$fileIdx == 1) + 3, nrow(cont))
    cont <- rhdf5::h5ls(res@h5files[2])
    expect_equal(sum(spdata$fileIdx == 2) + 3, nrow(cont))
})

test_that(".serialize_msfile_to_hdf5 works", {
    h5file <- tempfile()
    h5 <- rhdf5::H5Fcreate(h5file)
    rhdf5::h5createGroup(h5, "spectra")
    rhdf5::h5createGroup(h5, "checksum")
    rhdf5::H5Fclose(h5)
    md5 <- .serialize_msfile_to_hdf5(fileNames(sciex)[1], h5file)
    cont <- rhdf5::h5ls(h5file)
    expect_equal(nrow(cont), sum(fromFile(sciex) == 1) + 3)
})

test_that(".h5_read_bare works", {
    expect_error(.h5_read_bare())
    expect_error(.h5_read_bare("5"))

    fid <- .Call("_H5Fopen", sciex_h5@backend@h5files[1], 0L, PACKAGE = "rhdf5")
    res <- .h5_read_bare(fid, "/spectra/3")
    .Call("_H5Fclose", fid, PACKAGE = "rhdf5")
    res_2 <- rhdf5::h5read(sciex_h5@backend@h5files[1],
                           name = "/spectra/3")
    expect_equal(res, res_2)
    expect_equal(res[, 1], mz(sciex[[3]]))
    expect_equal(res[, 2], intensity(sciex[[3]]))
})

test_that(".h5_read_spectra works", {
    spd <- sciex_h5@spectraData
    res_h5 <- .h5_read_spectra(spd[spd$fileIdx == 1, , drop = FALSE],
                               sciex_h5@backend@h5files[1],
                               sciex_h5@backend@checksums[1])
    res_od <- spectra(sciex)[fromFile(sciex) == 1]
    expect_equal(res_h5, res_od)
    idx <- c(34, 65, 234, 453, 488)
    res_h5 <- .h5_read_spectra(spd[idx, ],
                               sciex_h5@backend@h5files[1],
                               sciex_h5@backend@checksums[1])
    expect_equal(res_h5, res_od[idx])
    ## MS1 & 2 data
    f <- msdata::proteomics(full.names = TRUE,
                            pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    tmt_h5 <- readMSnExperiment(f, path = tempdir(), backend = BackendHdf5())
    res_h5 <- .h5_read_spectra(tmt_h5@spectraData,
                               tmt_h5@backend@h5files,
                               tmt_h5@backend@checksums)
    res_od <- spectra(tmt_erwinia_on_disk)
    expect_equal(res_h5, res_od)
    expect_equal(res_od[123],
                 .h5_read_spectra(tmt_h5@spectraData[123, ],
                                  tmt_h5@backend@h5files,
                                  tmt_h5@backend@checksums))
})

test_that("backendReadSpectra,BackendHdf5 works", {
    sps <- spectra(sciex)
    be <- sciex_h5@backend
    spd <- sciex_h5@spectraData
    ## all spectra from the second file.
    res <- MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 2, ])
    expect_equal(res, sps[spd$fileIdx == 2])
    ## Some spectra
    idx <- c(23, 45, 12, 954, 976)
    res <- MSnbase:::backendReadSpectra(be, spd[idx, ])
    expect_equal(res, sps[idx])
})

test_that(".h5_write_spectra, and backendWriteSpectra,BackendHdf5 work", {
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
       system.file("microtofq/MM8.mzML", package = "msdata"))
    mse_h5 <- readMSnExperiment(f, path = paste0(tempdir(), "/1"),
                                backend = BackendHdf5())
    be <- mse_h5@backend
    sps <- spectra(mse_h5, return.type = "list")
    spd <- mse_h5@spectraData
    ## Write only 5 spectra to the second file.
    idx <- c(114:118)
    chksum <- MSnbase:::.h5_write_spectra(sps[idx], spd[idx, ], be@h5files[2])
    expect_true(chksum != be@checksums[2])
    res <- MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 1, ])
    expect_equal(sps[spd$fileIdx == 1], res)
    expect_error(MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 2, ]))
    expect_error(MSnbase:::backendReadSpectra(be, spd[spd$fileIdx == 2, ]))
    be@checksums[2] <- chksum
    res <- MSnbase:::backendReadSpectra(be, spd[idx, ])
    expect_equal(res, sps[idx])

    ## backendWriteSpectra, write spectra in arbitrary order
    idx <- c(which(spd$fileIdx == 1), idx)
    spd <- spd[idx, ]
    sps <- sps[idx]
    res <- MSnbase:::backendReadSpectra(be, spd)
    expect_equal(res, sps)

    idx <- c(34, 12, 5, 117, 114)
    chksum_orig <- be@checksums
    res <- MSnbase:::backendWriteSpectra(be, sps[idx], spd[idx, ])
    expect_true(is(res, "BackendHdf5"))
    expect_true(all(chksum_orig != res@checksums))
    res_sps <- MSnbase:::backendReadSpectra(res, spd[idx, ])
    expect_equal(res_sps, sps[idx])
})

test_that("backendSubset, BackendHdf5", {
    spd <- sciex_h5@spectraData
    spd <- spd[c(1000, 1003, 34, 64), ]
    res <- MSnbase:::backendSubset(sciex_h5@backend, spd)
    expect_equal(res@files, sciex_h5@backend@files[2:1])
    expect_equal(res@checksums, sciex_h5@backend@checksums[2:1])
    expect_equal(res@h5files, sciex_h5@backend@h5files[2:1])

    spd <- spd[3, , drop = FALSE]
    res <- MSnbase:::backendSubset(sciex_h5@backend, spd)
    expect_equal(res@files, sciex_h5@backend@files[1])
    expect_equal(res@checksums, sciex_h5@backend@checksums[1])
    expect_equal(res@h5files, sciex_h5@backend@h5files[1])
})
