test_that(".hdf5_group_name works", {
    res <- .hdf5_group_name(c("a/bb b/ccc", "a/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
    res <- .hdf5_group_name(
        c("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bb b/ccc",
          "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
})

sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
h5_sciex <- readHdf5DiskMSData(sf, hdf5file = tempfile())
f <- msdata::proteomics(
                 full.names = TRUE,
                 pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
h5_tmt <- readHdf5DiskMSData(f, hdf5file = tempfile())

test_that(".h5read_bare works", {
    grps <- MSnbase:::.hdf5_group_name(fileNames(h5_sciex))
    expect_error(.h5read_bare())
    expect_error(.h5read_bare("5"))
    res <- MSnbase:::.h5read_bare(h5_sciex@hdf5file, paste0(grps[1], "/F1.S003"))
    res_2 <- rhdf5::h5read(h5_sciex@hdf5file, paste0(grps[1], "/F1.S003"))
    expect_equal(res, res_2)
})

test_that(".read_spectra_hdf5 works", {
    res_h5 <- MSnbase:::.hdf5_read_spectra(fData(h5_sciex), h5_sciex@hdf5file,
                                           fileNames(h5_sciex))
    res_od <- spectra(sciex)
    expect_equal(res_h5, res_od)
    idx <- c(34, 65, 234, 453, 488)
    res_h5 <- MSnbase:::.hdf5_read_spectra(fData(h5_sciex)[idx, ], h5_sciex@hdf5file,
                                 fileNames(h5_sciex))
    expect_equal(res_h5, res_od[idx])
    ## MS1 & 2 data
    res_h5 <- MSnbase:::.hdf5_read_spectra(fData(h5_tmt), h5_tmt@hdf5file,
                                 fileNames(h5_tmt))
    res_od <- spectra(tmt_erwinia_on_disk)
    expect_equal(res_h5, res_od)
    expect_equal(res_od[123], MSnbase:::.hdf5_read_spectra(fData(h5_tmt)[123, ],
                                                 h5_tmt@hdf5file,
                                                 fileNames(h5_tmt)))
})

test_that(".apply_processing_queue works", {
    sps <- spectra(sciex)
    expect_equal(sps, MSnbase:::.apply_processing_queue(sps))
    q <- list(ProcessingStep("fromFile"))
    res <- MSnbase:::.apply_processing_queue(sps, q)
    expect_equal(res, lapply(sps, fromFile))

    q <- list(ProcessingStep(FUN = function(x, a) {
        fromFile(x) * a
    }, ARGS = list(a = 4)))
    res <- MSnbase:::.apply_processing_queue(sps, q)
    expect_equal(unlist(res), unlist(lapply(sps, fromFile)) * 4)

    q <- list(ProcessingStep("removePeaks", list(t = 0)),
              ProcessingStep("clean", list(all = TRUE)),
              ProcessingStep("intensity"))
    res <- MSnbase:::.apply_processing_queue(sps[1:3], q)
    expect_equal(res[[1]], intensity(clean(removePeaks(sps[[1]], t = 20),
                                           all = TRUE)))
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
