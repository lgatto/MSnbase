test_that("writeMSData works", {
    mzML_xsd_idx <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.2_idx.xsd",
                                                  package = "mzR"),
                                      isSchema = TRUE, useInternal = TRUE)

    ## using the onDisk data.
    odf <- tmt_erwinia_on_disk
    ## 1) Filter MS level 1, write, read and compare with tmt_erwinia_in_mem_ms1
    odf_out <- filterMsLevel(odf, msLevel = 1)
    out_file <- paste0(tempfile(), ".mzML")
    MSnbase:::.writeSingleMSData(odf_out, file = out_file,
                                 outformat = "mzml", copy = TRUE)
    ## Validating the mzML file
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    expect_equal(res$status, 0)
    odf_in <- readMSData(out_file, mode = "onDisk")
    ## Some stuff is different, i.e. totIonCurrent, basePeakMZ, basePeakIntensity
    expect_equal(unname(rtime(odf_in)), unname(rtime(tmt_erwinia_in_mem_ms1)))
    expect_equal(unname(mz(odf_in)), unname(mz(tmt_erwinia_in_mem_ms1)))
    expect_equal(unname(intensity(odf_in)),
                 unname(intensity(tmt_erwinia_in_mem_ms1)))
    ## feature data should be the same, expect spIdx, seqNum, spectrum and the
    ## columns re-calculated prior to saving.
    fd_out <- fData(odf_out)
    fd_in <- fData(odf_in)
    not_equal <- c("spIdx", "seqNum", "totIonCurrent", "basePeakMZ",
                   "basePeakIntensity", "spectrum")
    check_cols <- colnames(fd_out)[!(colnames(fd_out) %in% not_equal)]
    rownames(fd_out) <- NULL
    rownames(fd_in) <- NULL
    expect_equal(fd_in[, check_cols], fd_out[, check_cols])
    ## Explitely compare filterString
    expect_true(all(!is.na(fd_in$filterString)))
    expect_equal(fd_in$filterString, fd_out$filterString)
    ## 2) MS level 2, mzXML
    out_file <- paste0(tempfile(), ".mzXML")
    odf_out <- filterMsLevel(odf, msLevel = 2)
    expect_warning(
        MSnbase:::.writeSingleMSData(odf_out, file = out_file,
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
    ## Ensure that precursor data is preserved.
    fd_out <- fData(odf_out)
    fd_in <- fData(odf_in)
    rownames(fd_out) <- NULL
    rownames(fd_in) <- NULL
    ## Also some additional columns are different for mzXML output;
    ## filterString is also not supported by mzXML
    check_cols <- colnames(fd_out)[!(colnames(fd_out) %in%
                                     c(not_equal, "retentionTime",
                                       "precursorScanNum", "acquisitionNum",
                                       "injectionTime", "spectrumId",
                                       "filterString",
                                       "isolationWindowTargetMZ",
                                       "scanWindowLowerLimit",
                                       "scanWindowUpperLimit"))]
    expect_equal(fd_out[, check_cols], fd_in[, check_cols])
    ## Again force check:
    expect_equal(unname(precursorCharge(odf_out)),
                 unname(precursorCharge(odf_in)))
    expect_equal(unname(precursorMz(odf_out)),
                 unname(precursorMz(odf_in)))
    expect_equal(unname(precursorIntensity(odf_out)),
                 unname(precursorIntensity(odf_in)))

    ## Write two files.
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a.mzML", "/b.mzML"))
    MSnbase:::.writeMSData(microtofq_on_disk_ms1, file = out_file,
                           copy = FALSE)
    ## Validating the mzML file
    doc <- XML::xmlInternalTreeParse(out_file[1])
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    expect_equal(res$status, 0)
    doc <- XML::xmlInternalTreeParse(out_file[2])
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    expect_equal(res$status, 0)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))
    ## Providing software_processing.
    odf <- extdata_mzXML_on_disk
    out_file <- paste0(out_path, "/mzXML.mzML")
    MSnbase:::.writeMSData(odf, file = out_file,
                 software_processing = c("dummysoft", "0.0.1", "MS:1000035"))
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf)), unname(rtime(odf_in)))
    expect_equal(unname(mz(odf)), unname(mz(odf_in)))
    expect_equal(unname(intensity(odf)), unname(intensity(odf_in)))
})

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
    odf_proc <- filterMsLevel(tmt_erwinia_on_disk, msLevel = 1)
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

test_that("writeMSData,OnDiskMSnExp works", {
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a2.mzML", "/b2.mzML"))
    writeMSData(microtofq_on_disk_ms1, file = out_file, copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))

    ## Write MS1 and MS2
    out_file <- paste0(tempfile(), ".mzML")
    out_data <- tmt_erwinia_on_disk
    writeMSData(out_data, file = out_file, copy = FALSE)
    in_data <- readMSData(out_file, mode = "onDisk")
    expect_equal(rtime(out_data), rtime(in_data))
    ## Columns expected to be different:
    not_equal <- c("totIonCurrent", "basePeakMZ", "basePeakIntensity",
                   "lowMZ", "highMZ")
    check_cols <- !(colnames(fData(out_data)) %in% not_equal)
    ## without copying the data the actual numbers are different, but assignment
    ## is expected to be the same
    fData(out_data)$acquisitionNum <-
                      as.numeric(factor(fData(out_data)$acquisitionNum))
    fData(out_data)$precursorScanNum <-
                      as.numeric(factor(fData(out_data)$precursorScanNum))
    fData(in_data)$acquisitionNum <-
                      as.numeric(factor(fData(in_data)$acquisitionNum))
    fData(in_data)$precursorScanNum <-
                      as.numeric(factor(fData(in_data)$precursorScanNum))
    expect_equal(fData(out_data)[, check_cols], fData(in_data)[, check_cols])
    expect_equal(fData(out_data)$filterString, fData(in_data)$filterString)

    ## With copy = TRUE
    out_file <- paste0(tempfile(), ".mzML")
    out_data <- tmt_erwinia_on_disk
    writeMSData(out_data, file = out_file, copy = TRUE)
    in_data <- readMSData(out_file, mode = "onDisk")
    expect_equal(rtime(out_data), rtime(in_data))
    ## Columns expected to be different:
    not_equal <- c("totIonCurrent", "basePeakMZ", "basePeakIntensity",
                   "lowMZ", "highMZ")
    check_cols <- !(colnames(fData(out_data)) %in% not_equal)
    expect_equal(fData(out_data)[, check_cols], fData(in_data)[, check_cols])
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
    data_out <- readMSData(in_file, msLevel = 3, mode = "inMem")
    out_file <- paste0(tempfile(), ".mzML")
    writeMSData(data_out, file = out_file, outformat = "mzml", copy = TRUE)
    data_in <- readMSData(out_file, mode = "inMem", msLevel = 3)
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
    ## data_out <- readMSData(in_file, mode = "inMem", msLevel = 1)
    data_out <- as(data_out, "MSnExp")
    out_file <- paste0(tempfile(), ".mzML")
    writeMSData(data_out, file = out_file, outformat = "mzml", copy = FALSE)
    ## Reading the data as onDisk in, since we just compare the data anyway.
    data_in <- readMSData(out_file, mode = "onDisk", msLevel = 1)
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
