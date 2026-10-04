context("readMSData onDisk mode")

mzf <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
         MsDataHub::X20171016_POOL_POS_3_105.134.mzML())

## Load the data with readMSData mode = onDisk
odmse <- readMSData(files = mzf, centroided. = TRUE, mode = "onDisk")

test_that("Constructor performance and test for MSn", {
    ## Get the test data file.
    mzf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
    odmsn <- readMSData(files = mzf, centroided. = TRUE, mode = "onDisk")
    featDat <- fData(odmsn)
    featDat <- featDat[featDat$fileIdx == 1, ]
    ## Compare the constructors, i.e. the "new" and the C-level one for Spectrum1 and Spectrum2.
    ## system.time(
    ##     spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmsn))
    ## ) ## 4.1 sec.
    system.time(
        spM <- MSnbase:::.applyFun2SpectraOfFileMulti(featDat, filenames=fileNames(odmsn))
    ) ## 1.6 sec.
    ## expect_equal(spR, spM)
})
