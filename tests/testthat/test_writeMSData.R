test_that(".pattern_to_cv works", {
    ## Not found.
    expect_equal(.pattern_to_cv("unknown"), NA_character_)
    expect_equal(.pattern_to_cv("peak picking"), "MS:1000035")
    expect_equal(.pattern_to_cv("centroid"), "MS:1000035")
    expect_equal(.pattern_to_cv("Alignment/retention time adjustment"),
                 "MS:1000745")
})

test_that(".guessSoftwareProcessing works", {
    ## filterMsLevel: Filter: select MS level(s)
    odf_proc <- filterMsLevel(tmt_erwinia_on_disk, msLevel. = 1)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(res[[1]][1], "MSnbase")
    expect_equal(res[[1]][2], paste0(packageVersion("MSnbase"), collapse = "."))
    expect_equal(res[[1]][3], "MS:1002870")
    expect_equal(res[[1]][4], "MS:1001486")
    ## clean: Spectra cleaned NO CV YET
    ## bin: Spectra binned: NO CV YET
    ## removePeaks: Curves <= t set to '0': NO CV YET
    odf_proc <- removePeaks(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(res[[1]][4], "MS:1001486")
    expect_equal(length(res[[1]]), 4)
    ## normalise: Spectra normalised
    odf_proc <- normalise(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 5)
    expect_equal(res[[1]][4], "MS:1001486")
    expect_equal(res[[1]][5], "MS:1001484")
    ## pickPeaks: Spectra centroided
    odf_proc <- pickPeaks(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 6)
    expect_equal(res[[1]][6], "MS:1000035")
    ## smooth: Spectra smoothed
    odf_proc <- smooth(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 7)
    expect_equal(res[[1]][7], "MS:1000592")
    ## filterRt: Filter: select retention time
    odf_proc <- filterRt(odf_proc, rt = c(200, 600))
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 8)
    expect_equal(res[[1]][8], "MS:1001486")
    ## filterMz: FilteR: trim MZ
    ## filterFile: Filter: select file(s)
    ## filterAcquisitionNum: Filter: select by
    ## filterEmptySpectra: Removed XXX empty spectra
    ## And with providing in addition other processings.
    res <- .guessSoftwareProcessing(odf_proc, c("other_soft", "43.2.1"))
    expect_equal(res[[1]][8], "MS:1001486")
    expect_equal(res[[2]], c("other_soft", "43.2.1"))
    ## Check that we don't get unknown CV parameter.
    odf_proc <- clean(tmt_erwinia_on_disk)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(res[[1]][1], "MSnbase")
    expect_equal(res[[1]][2], paste0(packageVersion("MSnbase"), collapse = "."))
    expect_true(length(res[[1]]) == 3)
})


test_that("writeMSData,MSnExp works", {
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a3.mzML", "/b3.mzML"))
    writeMSData(microtofq_in_mem_ms1, file = out_file, copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))

    out_file <- paste0(out_path, c("/mzxml3.mzML"))
    writeMSData(extdata_mzXML_in_mem_ms2, file = out_file, copy = FALSE)
    odf_in <- readMSData(out_file, mode = "inMem")
    ## Check that main data is the same
    expect_equal(unname(rtime(odf_in)), unname(rtime(extdata_mzXML_in_mem_ms2)))
    expect_equal(unname(mz(odf_in)), unname(mz(extdata_mzXML_in_mem_ms2)))
    expect_equal(unname(intensity(odf_in)),
                 unname(intensity(extdata_mzXML_in_mem_ms2)))
    ## Check that header is the same
    expect_equal(unname(precursorMz(odf_in)),
                 unname(precursorMz(extdata_mzXML_in_mem_ms2)))
    expect_equal(unname(precursorCharge(odf_in)),
                 unname(precursorCharge(extdata_mzXML_in_mem_ms2)))
    expect_equal(unname(precursorIntensity(odf_in)),
                 unname(precursorIntensity(extdata_mzXML_in_mem_ms2)))

    in_file <- system.file(package = "msdata",
                           "proteomics/MS3TMT10_01022016_32917-33481.mzML.gz")
    data_out <- readMSData(in_file, msLevel. = 3, mode = "inMem")
    out_file <- paste0(tempfile(), ".mzML")
    writeMSData(data_out, file = out_file, outformat = "mzml", copy = TRUE)
    data_in <- readMSData(out_file, mode = "inMem", msLevel. = 3)
    expect_equal(rtime(data_in), rtime(data_out))
    expect_equal(mz(data_in), mz(data_out))
    expect_equal(intensity(data_in), intensity(data_out))
    expect_equal(precursorCharge(data_in), precursorCharge(data_out))
    expect_equal(precursorMz(data_in), precursorMz(data_out))
    expect_equal(precursorIntensity(data_in), precursorIntensity(data_out))
})

test_that("writeMSData works on CDF files", {
    in_file <- system.file(package = "msdata", "cdf/ko15.CDF")
    ## on disk
    data_out <- readMSData(in_file, mode = "onDisk")
    out_file <- paste0(tempfile(), ".mzML")
    writeMSData(data_out, file = out_file, outformat = "mzml", copy = FALSE)
    data_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(rtime(data_out), rtime(data_in))
    ## Most of the data in fData will be different, because CDF files do not
    ## provide that many header informations.
    sps_in <- spectra(data_in)
    sps_out <- spectra(data_out)
    expect_equal(lapply(sps_in, mz), lapply(sps_out, mz))
    expect_equal(lapply(sps_in, intensity), lapply(sps_out, intensity))

    ## in mem
    ## NOTE: reading a CDF file inMem is much slower than converting an
    ## OnDiskMSnExp into a MSnExp.
    ## data_out <- readMSData(in_file, mode = "inMem", msLevel. = 1)
    data_out <- as(data_out, "MSnExp")
    out_file <- paste0(tempfile(), ".mzML")
    writeMSData(data_out, file = out_file, outformat = "mzml", copy = FALSE)
    ## Reading the data as onDisk in, since we just compare the data anyway.
    data_in <- readMSData(out_file, mode = "onDisk", msLevel. = 1)
    expect_equal(rtime(data_out), rtime(data_in))
    expect_equal(mz(data_out), mz(data_in))
    expect_equal(intensity(data_out), intensity(data_in))
    sps_in <- spectra(data_in)
    sps_out <- spectra(data_out)
    tmp_fun <- function(z) {
        z@polarity <- 1L
        z
    }
    sps_in <- lapply(sps_in, tmp_fun)
    sps_out <- lapply(sps_out, tmp_fun)
    expect_equal(sps_in, sps_out)
})
