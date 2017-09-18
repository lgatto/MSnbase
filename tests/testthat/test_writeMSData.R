test_that("writeMSData works", {
    ## using the onDisk data.
    out_path <- tempdir()
    odf <- tmt_erwinia_on_disk
    ## filter MS level 1, write, read and compare with tmt_erwinia_in_mem_ms1
    odf_1 <- filterMsLevel(odf, msLevel = 1)
    out_file <- paste0(tempfile(), ".mzML")
    MSnbase:::.writeSingleMSData(odf_1, file = out_file,
                                 outformat = "mzml", copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    ## Some stuff is different, i.e. totIonCurrent, basePeakMZ, basePeakIntensity
    expect_equal(unname(rtime(odf_in)), unname(rtime(tmt_erwinia_in_mem_ms1)))
    expect_equal(unname(mz(odf_in)), unname(mz(tmt_erwinia_in_mem_ms1)))
    expect_equal(unname(intensity(odf_in)),
                 unname(intensity(tmt_erwinia_in_mem_ms1)))
    ## MS level 2, mzXML
    out_file <- paste0(out_path, "/write_test.mzXML")
    odf_2 <- filterMsLevel(odf, msLevel = 2)
    expect_warning(
        MSnbase:::.writeSingleMSData(odf_2, file = out_file,
                                     outformat = "mzxml", copy = TRUE)
    )
    odf_in <- readMSData(out_file, mode = "onDisk")
    ## retention time is saved with less precision in mzXML
    expect_true(all.equal(unname(rtime(odf_in)),
                          unname(rtime(tmt_erwinia_in_mem_ms2)),
                          tolerance = 1e4))
    expect_equal(unname(mz(odf_in)), unname(mz(tmt_erwinia_in_mem_ms2)))
    expect_equal(unname(intensity(odf_in)),
                 unname(intensity(tmt_erwinia_in_mem_ms2)))
    ## Write two files.
    out_file <- paste0(out_path, c("/a.mzML", "/b.mzML"))
    MSnbase:::.writeMSData(microtofq_on_disk_ms1, files = out_file,
                           copy = FALSE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))
    ## Providing software_processing.
    odf <- extdata_mzXML_on_disk
    out_file <- paste0(out_path, "/mzXML.mzML")
    MSnbase:::.writeMSData(x = odf, files = out_file,
                 software_processing = c("dummysoft", "0.0.1", "MS:1000035"))
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf)), unname(rtime(odf_in)))
    expect_equal(unname(mz(odf)), unname(mz(odf_in)))
    expect_equal(unname(intensity(odf)), unname(intensity(odf_in)))
})

test_that(".pattern_to_cv works", {
    ## Not found.
    expect_equal(.pattern_to_cv("unknown"), "MS:-1")
    expect_equal(.pattern_to_cv("peak picking"), "MS:1000035")
    expect_equal(.pattern_to_cv("centroid"), "MS:1000035")    
    expect_equal(.pattern_to_cv("Alignment/retention time adjustment"),
                 "MS:1000745")
})

test_that(".guessSoftwareProcessing works", {
    ## filterMsLevel: Filter: select MS level(s)
    odf_proc <- filterMsLevel(tmt_erwinia_on_disk, msLevel = 1)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(res[[1]][1], "MSnbase")
    expect_equal(res[[1]][2], paste0(packageVersion("MSnbase"), collapse = "."))
    expect_equal(res[[1]][3], "MS:1001486")
    ## clean: Spectra cleaned NO CV YET
    ## bin: Spectra binned: NO CV YET
    ## removePeaks: Curves <= t set to '0': NO CV YET
    odf_proc <- removePeaks(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(res[[1]][3], "MS:1001486")
    expect_equal(length(res[[1]]), 3)
    ## normalise: Spectra normalised
    odf_proc <- normalise(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 4)
    expect_equal(res[[1]][3], "MS:1001486")
    expect_equal(res[[1]][4], "MS:1001484")
    ## pickPeaks: Spectra centroided
    odf_proc <- pickPeaks(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 5)
    expect_equal(res[[1]][5], "MS:1000035")
    ## smooth: Spectra smoothed
    odf_proc <- smooth(odf_proc)
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 6)
    expect_equal(res[[1]][6], "MS:1000542")
    ## filterRt: Filter: select retention time
    odf_proc <- filterRt(odf_proc, rt = c(200, 600))
    res <- .guessSoftwareProcessing(odf_proc)
    expect_equal(length(res[[1]]), 7)
    expect_equal(res[[1]][7], "MS:1001486")
    ## filterMz: FilteR: trim MZ
    ## filterFile: Filter: select file(s)
    ## filterAcquisitionNum: Filter: select by
    ## filterEmptySpectra: Removed XXX empty spectra
    ## And with providing in addition other processings.
    res <- .guessSoftwareProcessing(odf_proc, c("other_soft", "43.2.1"))
    expect_equal(res[[1]][7], "MS:1001486")
    expect_equal(res[[2]], c("other_soft", "43.2.1"))
})


test_that("write,OnDiskMSnExp works", {
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a2.mzML", "/b2.mzML"))
    write(microtofq_on_disk_ms1, files = out_file, copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))
})

test_that("write,MSnExp works", {
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a3.mzML", "/b3.mzML"))
    write(microtofq_in_mem_ms1, files = out_file, copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))

    out_file <- paste0(out_path, c("/mzxml3.mzML"))
    write(extdata_mzXML_in_mem_ms2, files = out_file, copy = FALSE)
    odf_in <- readMSData(out_file, mode = "inMem")
    expect_equal(unname(rtime(odf_in)), unname(rtime(extdata_mzXML_in_mem_ms2)))
    ## precursor scans are missing, thus we don't have precursor data in the
    ## output file.
    expect_true(all(precursorMz(odf_in) != precursorMz(extdata_mzXML_in_mem_ms2)))

    expect_equal(mz(odf_in), mz(extdata_mzXML_in_mem_ms2))
    expect_equal(intensity(odf_in), intensity(extdata_mzXML_in_mem_ms2))    
})
