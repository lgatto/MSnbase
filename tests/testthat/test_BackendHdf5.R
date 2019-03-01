context("BackendHdf5 class")

test_that("BackendHdf5 validators work", {
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
    tst@modCount <- 1L
    expect_error(validObject(tst))
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
                               sciex_h5@backend@modCount[1])
    res_od <- spectra(sciex)[fromFile(sciex) == 1]
    expect_equal(res_h5, res_od)
    idx <- c(34, 65, 234, 453, 488)
    res_h5 <- .h5_read_spectra(spd[idx, ],
                               sciex_h5@backend@h5files[1],
                               sciex_h5@backend@modCount[1])
    expect_equal(res_h5, res_od[idx])
    ## MS1 & 2 data
    f <- msdata::proteomics(full.names = TRUE,
                            pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    tmt_h5 <- readMSnExperiment(f, path = tempdir(), backend = BackendHdf5())
    res_h5 <- .h5_read_spectra(tmt_h5@spectraData,
                               tmt_h5@backend@h5files,
                               tmt_h5@backend@modCount)
    res_od <- spectra(tmt_erwinia_on_disk)
    expect_equal(res_h5, res_od)
    expect_equal(res_od[123],
                 .h5_read_spectra(tmt_h5@spectraData[123, ],
                                  tmt_h5@backend@h5files,
                                  tmt_h5@backend@modCount))
})

test_that("backendReadSpectra,BackendHdf5 works", {
    sps <- spectra(sciex)
    be <- sciex_h5@backend
    spd <- sciex_h5@spectraData
    ## all spectra from the second file.
    res <- backendReadSpectra(be, spd[spd$fileIdx == 2, ])
    expect_equal(res, sps[spd$fileIdx == 2])
    ## Some spectra
    idx <- c(23, 45, 12, 954, 976)
    res <- backendReadSpectra(be, spd[idx, ])
    expect_equal(res, sps[idx])
})

test_that(".h5_write_spectra, and backendWriteSpectra,BackendHdf5 work", {
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
       system.file("microtofq/MM8.mzML", package = "msdata"))
    mse_h5 <- readMSnExperiment(f, path = paste0(tempdir(), "/1"),
                                backend = BackendHdf5())
    be <- mse_h5@backend
    sps <- spectrapply(mse_h5)
    spd <- mse_h5@spectraData
    ## Write only 5 spectra to the second file.
    idx <- c(114:118)
    .h5_write_spectra(sps[idx], spd[idx, ], be@h5files[2], be@modCount[2] + 1L)
    res <- backendReadSpectra(be, spd[spd$fileIdx == 1, ])
    expect_equal(sps[spd$fileIdx == 1], res)
    expect_error(backendReadSpectra(be, spd[spd$fileIdx == 2, ]),
                 "The data .* have changed")
    be@modCount[2] <- 1L
    res <- backendReadSpectra(be, spd[idx, ])
    expect_equal(res, sps[idx])

    ## backendWriteSpectra, write spectra in arbitrary order
    idx <- c(which(spd$fileIdx == 1), idx)
    spd <- spd[idx, ]
    sps <- sps[idx]
    res <- backendReadSpectra(be, spd)
    expect_equal(res, sps)

    idx <- c(34, 12, 5, 117, 114)
    modCount_orig <- be@modCount
    res <- backendWriteSpectra(be, sps[idx], spd[idx, ])
    expect_true(is(res, "BackendHdf5"))
    expect_true(all(modCount_orig != res@modCount))
    res_sps <- backendReadSpectra(res, spd[idx, ])
    expect_equal(res_sps, sps[idx])
})

test_that("backendSubset,BackendHdf5", {
    spd <- sciex_h5@spectraData
    spd <- spd[c(1000, 1003, 34, 64), ]
    res <- backendSubset(sciex_h5@backend, spd)
    expect_equal(res@files, sciex_h5@backend@files[2:1])
    expect_equal(res@modCount, sciex_h5@backend@modCount[2:1])
    expect_equal(res@h5files, sciex_h5@backend@h5files[2:1])

    spd <- spd[3, , drop = FALSE]
    res <- backendSubset(sciex_h5@backend, spd)
    expect_equal(res@files, sciex_h5@backend@files[1])
    expect_equal(res@modCount, sciex_h5@backend@modCount[1])
    expect_equal(res@h5files, sciex_h5@backend@h5files[1])
})

test_that("backendSplitByFile,BackendHdf5 works", {
    b <- BackendHdf5()
    b@files <- c("a", "b", "c")
    tmpfiles <- paste0(tempfile(), letters[1:4], ".h5")
    file.create(tmpfiles)
    on.exit(unlink(tmpfiles))
    b@h5files <- tmpfiles[1:3]
    b@modCount <- rep(0L, 3L)
    spd <- DataFrame(fileIdx = c(3, 3, 1, 1, 2, 1))
    res <- backendSplitByFile(b, spd)
    bl <- BackendHdf5()
    bl@files <- "a"
    bl@h5files <- tmpfiles[1]
    bl@modCount <- 0L
    l <- list("1"=bl, "2"=bl, "3"=bl)
    l[[2]]@files <- "b"
    l[[2]]@h5files <- tmpfiles[2]
    l[[3]]@files <- "c"
    l[[3]]@h5files <- tmpfiles[3]
    expect_equal(backendSplitByFile(b, spd), l)
    r <- b
    r@files[1] <- "d"
    r@h5files[1] <- tmpfiles[4]
    r@modCount[1L] <- 1L
    l[[1]]@files <- "d"
    l[[1]]@h5files <- tmpfiles[4]
    l[[1]]@modCount <- 1L
    backendSplitByFile(b, spd) <- l
    expect_equal(b, r)
})
