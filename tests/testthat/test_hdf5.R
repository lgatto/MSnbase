sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
h5_sciex <- readHdf5DiskMSData(sf, hdf5path = tempdir())
f <- msdata::proteomics(
                 full.names = TRUE,
                 pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
h5_tmt <- readHdf5DiskMSData(f, hdf5path = tempdir())

test_that("validHdf5MSnExp works", {
    expect_true(MSnbase:::validHdf5MSnExp(h5_sciex))
    expect_true(validObject(h5_sciex))
    tmp <- h5_sciex
    tmp@hdf5file <- tmp@hdf5file[2]
    res <- MSnbase:::validHdf5MSnExp(tmp, check_md5 = FALSE)
    expect_true(is.character(res))
    tmp@hdf5file <- c(tmp@hdf5file, "a")
    res <- MSnbase:::validHdf5MSnExp(tmp, check_md5 = FALSE)
    expect_true(is.character(res))
    expect_error(validObject(tmp))

    expect_true(MSnbase:::.h5_check_md5(h5_sciex))
})

test_that(".serialize_msfile_to_hdf5 works", {
    h5 <- tempfile()
    md5 <- MSnbase:::.serialize_msfile_to_hdf5(fileNames(sciex)[1], h5)
    cont <- rhdf5::h5ls(h5)
    expect_equal(nrow(cont), sum(fromFile(sciex) == 1) + 2)
})

test_that("serialize_to_hdf5 works", {
})

test_that("readHdf5DiskMSnData works", {
})

test_that("hdf5FileName works", {
})

test_that(".h5_group_name works", {
    res <- MSnbase:::.h5_group_name(c("a/bb b/ccc", "a/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
    res <- MSnbase:::.h5_group_name(
        c("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bb b/ccc",
          "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
})


test_that(".h5_read_bare works", {
    grps <- MSnbase:::.h5_group_name(fileNames(h5_sciex))
    expect_error(.h5_read_bare())
    expect_error(.h5_read_bare("5"))
    fid <- .Call("_H5Fopen", h5_sciex@hdf5file, 0L, PACKAGE = "rhdf5")
    res <- MSnbase:::.h5_read_bare(fid, paste0(grps[1], "/3"))
    .Call("_H5Fclose", fid, PACKAGE = "rhdf5")
    res_2 <- rhdf5::h5read(h5_sciex@hdf5file[1],
                           name = paste0(grps[1], "/3"))
    expect_equal(res, res_2)
    expect_equal(res[, 1], mz(sciex[[3]]))
    expect_equal(res[, 2], intensity(sciex[[3]]))
})

test_that(".read_spectra_hdf5 works", {
    fd <- fData(h5_sciex)
    res_h5 <- MSnbase:::.h5_read_spectra(fd[fd$fileIdx == 1, , drop = FALSE],
                                         h5_sciex@hdf5file[1],
                                         h5_sciex@md5sum[1],
                                         fileNames(h5_sciex)[1])
    res_od <- spectra(sciex)[fromFile(sciex) == 1]
    expect_equal(res_h5, res_od)
    idx <- c(34, 65, 234, 453, 488)
    res_h5 <- MSnbase:::.h5_read_spectra(fd[idx, ], h5_sciex@hdf5file[1],
                                         h5_sciex@md5sum[1],
                                         fileNames(h5_sciex)[1])
    expect_equal(res_h5, res_od[idx])
    ## MS1 & 2 data
    res_h5 <- MSnbase:::.h5_read_spectra(fData(h5_tmt), h5_tmt@hdf5file,
                                         h5_tmt@md5sum,
                                         fileNames(h5_tmt))
    res_od <- spectra(tmt_erwinia_on_disk)
    expect_equal(res_h5, res_od)
    expect_equal(res_od[123], MSnbase:::.h5_read_spectra(fData(h5_tmt)[123, ],
                                                         h5_tmt@hdf5file,
                                                         h5_tmt@md5sum,
                                                         fileNames(h5_tmt)))
})

test_that("spectrapply,Hdf5MSnExp works", {
    library(BiocParallel)
    sps <- spectra(sciex)
    h5_sps <- spectrapply(h5_sciex)
    expect_equal(sps, h5_sps)

    register(SerialParam())
    h5_sps <- spectrapply(h5_sciex, BPPARAM = SerialParam())
    expect_equal(sps, h5_sps)

    tmp_tmt <- smooth(tmt_erwinia_on_disk, halfWindowSize = 4L)
    tmp_h5 <- smooth(h5_tmt, halfWindowSize = 4L)
    expect_warning(expect_equal(spectra(tmp_tmt), spectra(tmp_h5)))
})

test_that(".h5_write_spectra works", {
    one_f <- filterFile(sciex, 1)
    sps <- spectra(one_f)
    fd <- fData(one_f)
    names(sps) <- fd$spIdx
    h5f <- tempfile()
    h5 <- H5Fcreate(h5f)
    H5Fclose(h5)
    fln <- "test"
    grp <- MSnbase:::.h5_group_name(fln)
    md5 <- MSnbase:::.h5_write_spectra(sps, h5f, group = grp, prune = FALSE)

    expect_error(MSnbase:::.h5_read_spectra(fd, h5f, "123", fln))
    res <- MSnbase:::.h5_read_spectra(fd, h5f, md5, fln)
    expect_equal(unname(res), unname(sps))

    ## Update writing only a subset.
    idx <- c(1, 45, 65, 123, 434, 894)
    md5 <- MSnbase:::.h5_write_spectra(sps[idx], h5f, group = grp, prune = TRUE)
    res <- MSnbase:::.h5_read_spectra(fd[idx, ], h5f, md5, fln)
    expect_equal(unname(res), unname(sps[idx]))
    cont <- rhdf5::h5ls(h5f)
    expect_equal(nrow(cont), length(idx) + 2)

    ## Writing empty spectra?
    sps[[1]] <- clean(removePeaks(sps[[1]], t = 10e9), all = TRUE)
    idx <- c(1, 12, 13, 14)
    md5 <- MSnbase:::.h5_write_spectra(sps[idx], h5f, group = grp, prune = TRUE)
    res <- MSnbase:::.h5_read_spectra(fd[idx, ], h5f, md5, fln)
    res[[1]]@tic <- 0 # fData still contains the TIC, while the spectra has 0.
    expect_equal(unname(res), unname(sps[idx]))
    cont <- rhdf5::h5ls(h5f)
    expect_equal(nrow(cont), length(idx) + 2)
})

test_that("filterFile,Hdf5MSnExp works", {
    res <- filterFile(h5_sciex)
    expect_equal(fileNames(res), fileNames(h5_sciex))
    expect_equal(res@hdf5file, h5_sciex@hdf5file)
    expect_equal(fData(res), fData(h5_sciex))

    res <- filterFile(h5_sciex, 1)
    expect_equal(fileNames(res), fileNames(h5_sciex)[1])
    expect_equal(res@hdf5file, h5_sciex@hdf5file[1])
    fd <- fData(h5_sciex)
    expect_equal(fData(res), fd[fd$fileIdx == 1, ])

    res <- filterFile(h5_sciex, 2)
    expect_equal(fileNames(res), fileNames(h5_sciex)[2])
    expect_equal(res@hdf5file, h5_sciex@hdf5file[2])
    fd <- fData(h5_sciex)
    fd <- fd[fd$fileIdx == 2, ]
    fd$fileIdx <- 1
    expect_equal(fData(res), fd)
})

test_that("writeHdf5Data works", {
    ## Filter by retention time.
    sciex_flt <- filterRt(h5_sciex, rt = c(10, 100))
    sps <- spectra(sciex_flt)
    ## consolidate will update the hdf5 files, so, h5_sciex will fail!
    res <- MSnbase:::.consolidate(sciex_flt)
    ##
    expect_true(MSnbase:::.h5_check_md5(res))
    expect_false(MSnbase:::.h5_check_md5(h5_sciex))
    expect_error(spectra(h5_sciex))
    cont <- rhdf5::h5ls(res@hdf5file[1])
    expect_equal(nrow(cont), (length(sps) / 2) + 2)

    ## Doing some data manipulations on the subsetted data
    res_cent <- pickPeaks(smooth(res))
    expect_equal(length(res_cent@spectraProcessingQueue), 2)
    sps_cent <- spectra(res_cent)
    sps_cent <- lapply(sps_cent, function(z) {
        z@tic <- 0
        z
    })
    ## Consolidate with in-place replacement.
    md5_before <- res_cent@md5sum
    writeHdf5Data(res_cent)
    md5_after <- res_cent@md5sum
    expect_true(all(md5_before != md5_after))

    # Don't have a processing queue because we replaced in-place
    expect_equal(length(res_cent@spectraProcessingQueue), 0)
    res_cent_sps <- spectra(res_cent)
    res_cent_sps <- lapply(res_cent_sps, function(z) {
        z@tic <- 0
        z
    })
    ## The spectra after consilidation are the same as the one before
    ## consolidation that used the on-the-fly data manipulations
    expect_equal(sps_cent, res_cent_sps)
    expect_true(peaksCount(res_cent_sps[[1]]) < peaksCount(sps[[1]]))

    ## Restore the hdf5 files again.
    file.remove(h5_sciex@hdf5file)
    sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    h5_sciex <- readHdf5DiskMSData(sf, hdf5path = tempdir())
})
