test_that("writeMSData works", {
    mzML_xsd_idx <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.2_idx.xsd",
                                                  package = "mzR"),
                                      isSchema = TRUE, useInternalNodes = TRUE)

    ## using the onDisk data.
    odf <- tmt_erwinia_on_disk
    ## 1) Filter MS level 1, write, read and compare with tmt_erwinia_in_mem_ms1
    odf_out <- filterMsLevel(odf, msLevel. = 1)
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
    odf_out <- filterMsLevel(odf, msLevel. = 2)
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



test_that("writeMSData,OnDiskMSnExp works", {
    out_path <- tempdir()
    out_file <- paste0(out_path, c("/a2.mzML", "/b2.mzML"))
    writeMSData(microtofq_on_disk_ms1, file = out_file, copy = TRUE)
    odf_in <- readMSData(out_file, mode = "onDisk")
    expect_equal(unname(rtime(odf_in)), unname(rtime(microtofq_in_mem_ms1)))
    expect_equal(spectra(odf_in), spectra(microtofq_in_mem_ms1))
    ##
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
    ##
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