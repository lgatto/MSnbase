test_that("msLevel set correctly", {
    f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
    x <- readMSData(f, verbose = FALSE, mode = "onDisk")
    ms1 <- filterMsLevel(x, msLevel = 1)
    ms2 <- filterMsLevel(x, msLevel = 2)
    expect_true(all(!centroided(ms1)))
    expect_true(all(centroided(ms2)))
    x <- readMSData(f, centroided. = TRUE, verbose = FALSE, mode = "onDisk")
    expect_true(all(centroided(x)))
    x <- readMSData(f, centroided. = FALSE, verbose = FALSE, mode = "onDisk")
    expect_true(all(!centroided(x)))
    x <- readMSData(f, centroided. = c(FALSE, TRUE), verbose = FALSE, mode = "onDisk")
    expect_true(all(centroided(filterMsLevel(x, msLevel = 2))))
    expect_true(all(!centroided(filterMsLevel(x, msLevel = 1))))
    x2 <- readMSData(f, centroided. = c(FALSE, TRUE, NA, NA), verbose = FALSE, mode = "onDisk")
    expect_identical(centroided(x), centroided(x2))
    ## In mem with centroided.
    f <- system.file("microtofq/MM14.mzML", package = "msdata")
    x <- readMSData(f, msLevel = 1)
    expect_true(all(centroided(x)))
    x <- readMSData(f, msLevel. = 1, centroided. = FALSE)
    expect_true(all(!centroided(x)))
})

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

test_that("readMSData inMemory and onDisk reading CDF", {
    library(msdata)
    f <- system.file("cdf/ko15.CDF",  package = "msdata")
    odmse <- readMSData(f, mode = "onDisk")
    mse <- readMSData(f, msLevel. = 1, mode = "inMemory")
    all.equal(spectra(odmse), spectra(mse))
})
