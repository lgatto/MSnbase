context("readMSData2")

############################################################
## Load the required data files.
.getMzMLFiles <- function(force.msdata = FALSE) {
    ## Return the mzML files, the ones from the XXX package, or if run
    ## locally, some of my test files.
    HOST <- unlist(strsplit(system("hostname", intern = TRUE), split = ".",
                            perl = FALSE, fixed = TRUE))[1]
    if (HOST == "macbookjo" & !force.msdata) {
        mzfiles <- dir("/Users/jo/R-workspaces/EURAC/2016/2016-04-21-PolarMetabolom/data/mzML/",
                       pattern = "POS_C_O", full.names = TRUE)
    } else {
        require(msdata)
        mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                     system.file("microtofq/MM8.mzML", package = "msdata"))
    }
    return(mzfiles)
}

mzf <- .getMzMLFiles()[1:2]

## Load the data with readMSData2
odmse <- MSnbase:::readMSData2(files = mzf, centroided = TRUE, backend = "disk")

test_that("Constructor performance and test for MS1 only", {
    featDat <- fData(odmse)
    featDat <- featDat[featDat$fileIdx == 1, ]
    system.time(
        spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmse))
    ) ## 19.5 sec.
    system.time(
        spM <- MSnbase:::.applyFun2SpectraOfFileMulti(featDat, filenames=fileNames(odmse))
    ) ## 3.2 sec.
    expect_equal(spR, spM)
})

test_that("Constructor performance and test for MSn", {
    ## Get the test data file.
    mzf <- system.file("proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                       package = "msdata")
    odmsn <- MSnbase:::readMSData2(files = mzf, centroided = TRUE)
    featDat <- fData(odmsn)
    featDat <- featDat[featDat$fileIdx == 1, ]
    ## Compare the constructors, i.e. the "new" and the C-level one for Spectrum1 and Spectrum2.
    system.time(
        spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmsn))
    ) ## 4.1 sec.
    system.time(
        spM <- MSnbase:::.applyFun2SpectraOfFileMulti(featDat, filenames=fileNames(odmsn))
    ) ## 1.6 sec.
    expect_equal(spR, spM)
})
