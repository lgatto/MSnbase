context("OnDiskMSnExp class, 2")

f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()

test_that("OnDiskMSnExp constructor", {
    x <- tmt_erwinia_on_disk
    expect_true(validObject(x))
    expect_true(all(unique(msLevel(x)) == 1:2))
    expect_true(all(isCurrent(x)))
    expect_true(isVersioned(x))
    expect_null(show(x))
    x1 <- tmt_erwinia_on_disk_ms1
    expect_true(validObject(x1))
    expect_true(unique(msLevel(x1)) == 1)
    x2 <- tmt_erwinia_on_disk_ms2
    expect_true(validObject(x2))
    expect_true(unique(msLevel(x2)) == 2)
    expect_identical(length(x), length(x1) + length(x2))
    expect_identical(featureNames(x),
                     sort(c(featureNames(x1), featureNames(x2))))
    expect_identical(fileNames(x), fileNames(x1))
    expect_identical(fileNames(x), fileNames(x2))
})



test_that("Write mgf", {
    x1 <- tmt_erwinia_in_mem_ms2
    x2 <- tmt_erwinia_on_disk_ms2
    tf1 <- tempfile()
    tf2 <- tempfile()
    writeMgfData(x1[2:4], con = tf1)
    writeMgfData(x2[2:4], con = tf2)
    mgf1 <- readMgfData(tf1)
    mgf2 <- readMgfData(tf2)
    mgf2@processingData <- mgf1@processingData
    mgf2@phenoData <- mgf1@phenoData
    expect_equal(mgf1, mgf2)
    unlink(tf1)
    unlink(tf2)
})

test_that("Adding identification data", {
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
    x1 <- extdata_mzXML_in_mem_ms2
    x1 <- addIdentificationData(x1, identFile, verbose = FALSE)
    x2 <- extdata_mzXML_on_disk_ms2
    x2 <- addIdentificationData(x2, identFile, verbose = FALSE)
    fv <- intersect(fvarLabels(x1), fvarLabels(x2))
    expect_identical(fData(x1)[, fv], fData(x2)[, fv])
    expect_identical(idSummary(x1), idSummary(x2))
    x1rm <- removeMultipleAssignment(x1)
    x2rm <- removeMultipleAssignment(x2)
    expect_identical(fData(x1rm)[, fv], fData(x2rm)[, fv])
    x1nid <- removeNoId(x1)
    x2nid <- removeNoId(x2)
    expect_identical(fData(x1nid)[, fv], fData(x2nid)[, fv])
})


test_that("Spectrum processing", {
    ## bin
    ## clean -> test_OnDiskMSnExp.R
    ## normalise/normalize -> test_OnDiskMSnExp_other_methods.R
    ## pickPeaks -> test_OnDiskMSnExp_other_methods.R
    ## quantify
    ## removePeaks -> test_OnDiskMSnExp.R
    ## removeReporters -> here
    ## smooth -> test_OnDiskMSnExp_other_methods.R
    ## trimMz -> test_OnDiskMSnExp_other_methods.R
})
