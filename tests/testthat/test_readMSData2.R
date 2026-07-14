context("readMSData onDisk mode")

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

mzf <- .getMzMLFiles(TRUE)[1:2]

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
