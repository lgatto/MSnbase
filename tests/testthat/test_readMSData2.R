context("readMSData onDisk mode")

test_that("msLevel set correctly", {
    f <- msdata::proteomics(full.names = TRUE,
                            pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    x <- readMSData(f, verbose = FALSE, mode = "onDisk")
    ms1 <- filterMsLevel(x, msLevel = 1)
    ms2 <- filterMsLevel(x, msLevel = 2)
    expect_true(all(!centroided(ms1)))
    expect_true(all(centroided(ms2)))
    x <- readMSData(f, centroided = TRUE, verbose = FALSE, mode = "onDisk")
    expect_true(all(centroided(x)))
    x <- readMSData(f, centroided = FALSE, verbose = FALSE, mode = "onDisk")
    expect_true(all(!centroided(x)))
    x <- readMSData(f, centroided = c(FALSE, TRUE), verbose = FALSE, mode = "onDisk")
    expect_true(all(centroided(filterMsLevel(x, msLevel = 2))))
    expect_true(all(!centroided(filterMsLevel(x, msLevel = 1))))    
    x2 <- readMSData(f, centroided = c(FALSE, TRUE, NA, NA), verbose = FALSE, mode = "onDisk")
    expect_identical(centroided(x), centroided(x2))
    ## In mem with centroided.
    f <- system.file("microtofq/MM14.mzML", package = "msdata")
    x <- readMSData(f, msLevel = 1)
    expect_true(all(centroided(x)))
    x <- readMSData(f, msLevel = 1, centroided = FALSE)
    expect_true(all(!centroided(x)))
})

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
odmse <- readMSData(files = mzf, centroided = TRUE, mode = "onDisk")

test_that("Constructor performance and test for MS1 only", {
    featDat <- fData(odmse)
    featDat <- featDat[featDat$fileIdx == 1, ]
    ## system.time(
    ##     spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmse))
    ## ) ## 19.5 sec.
    system.time(
        spM <- MSnbase:::.applyFun2SpectraOfFileMulti(featDat, filenames=fileNames(odmse))
    ) ## 3.2 sec.
    ## expect_equal(spR, spM)
})

test_that("Constructor performance and test for MSn", {
    ## Get the test data file.
    mzf <- proteomics(full.name = TRUE,
                      pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    odmsn <- readMSData(files = mzf, centroided = TRUE, mode = "onDisk")
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

test_that("readMSData inMemory and onDisk reading CDF", {
    library(msdata)
    f <- system.file("cdf/ko15.CDF",  package = "msdata")
    odmse <- readMSData(f, mode = "onDisk")
    mse <- readMSData(f, msLevel. = 1, mode = "inMemory")
    all.equal(spectra(odmse), spectra(mse))
})
