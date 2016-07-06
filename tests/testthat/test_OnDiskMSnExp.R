context("OnDiskMSnExp class")

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
## Load the data as an MSnExp into memory.
suppressWarnings(
    inMem <- readMSData(files = mzf, msLevel = 1, centroided = TRUE,
                        verbose = FALSE)
)
## Load the data as OnDiskMSnExp.
suppressWarnings(
    onDisk <- readMSData2(files = mzf, msLevel = 1, centroided=TRUE,
                          verbose = FALSE)
)

############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
    expect_identical(as.character(class(inMem)), "MSnExp")
    expect_identical(as.character(class(onDisk)), "OnDiskMSnExp")
    expect_true(validObject(onDisk))
})

############################################################
## compare MSnExp against OnDiskMSnExp
test_that("Compare MS1 MSnExp and OnDiskMSnExp content", {
    ## o Compare spectra values.
    expect_true(all.equal(inMem, onDisk))
    ## o fromFile
    expect_identical(fromFile(inMem), fromFile(onDisk))
    ## o msLevel
    expect_identical(msLevel(inMem), msLevel(onDisk))
    ## o acquisitionNum
    expect_identical(acquisitionNum(inMem), acquisitionNum(onDisk))
    ## o scanIndex
    expect_identical(scanIndex(inMem), scanIndex(onDisk))
    ## o centroided
    expect_identical(centroided(inMem), centroided(onDisk))
    centroided(inMem) <- FALSE
    centroided(onDisk) <- FALSE
    expect_identical(centroided(inMem), centroided(onDisk))
    expect_that(centroided(onDisk) <- c(TRUE, FALSE, TRUE), throws_error())
    ## o rtime
    expect_identical(rtime(inMem), rtime(onDisk))
    ## o polarity
    expect_identical(polarity(inMem), polarity(onDisk))
    ## o tic
    expect_identical(tic(inMem), tic(onDisk))
    ## o ionCount
    expect_identical(ionCount(inMem), ionCount(onDisk))
    ## o peaksCount
    expect_identical(peaksCount(inMem), peaksCount(onDisk))
    ## o intensity
    expect_identical(intensity(inMem), intensity(onDisk))
    ## o mz
    expect_identical(mz(inMem), mz(onDisk))
})

############################################################
## Compare cleaned data.
## o spectra.
## o ionCount.
## o tic
## o peaksCount.
test_that("Compare removePeaks and cleaned MSnExp and OnDiskMSnExp", {
    ## o clean
    inMemCleaned <- clean(inMem)
    onDiskCleaned <- clean(onDisk)
    expect_true(all.equal(inMemCleaned, onDiskCleaned))
    expect_identical(ionCount(inMemCleaned), ionCount(onDiskCleaned))
    expect_identical(tic(inMemCleaned), tic(onDiskCleaned))
    expect_identical(peaksCount(inMemCleaned), peaksCount(onDiskCleaned))

    ## o removePeaks
    inMemRemPeaks <- removePeaks(inMem, t = 1000)
    onDiskRemPeaks <- removePeaks(onDisk, t = 1000)
    expect_true(all.equal(inMemRemPeaks, onDiskRemPeaks))
    expect_identical(ionCount(inMemRemPeaks), ionCount(onDiskRemPeaks))
    expect_identical(tic(inMemRemPeaks), tic(onDiskRemPeaks))
    expect_identical(peaksCount(inMemRemPeaks), peaksCount(onDiskRemPeaks))

    ## o removePeaks and clean
    inMemRemPeaksCleaned <- clean(inMemRemPeaks)
    onDiskRemPeaksCleaned <- clean(onDiskRemPeaks)
    expect_true(all.equal(inMemRemPeaksCleaned, onDiskRemPeaksCleaned))
    expect_identical(ionCount(inMemRemPeaksCleaned),
                     ionCount(onDiskRemPeaksCleaned))
    expect_identical(tic(inMemRemPeaksCleaned), tic(onDiskRemPeaksCleaned))
    expect_identical(peaksCount(inMemRemPeaksCleaned),
                     peaksCount(onDiskRemPeaksCleaned))

    ## compare assayData, intensity and mz,
    expect_equal(assayData(inMemRemPeaksCleaned),
                 assayData(onDiskRemPeaksCleaned))
    expect_equal(intensity(inMemRemPeaksCleaned),
                 intensity(onDiskRemPeaksCleaned))
    expect_equal(mz(inMemRemPeaksCleaned), mz(onDiskRemPeaksCleaned))
})


############################################################
## [[
test_that("Compare subsetting between OnDiskMSnExp and MSnExp", {
    sp1 <- inMem[[77]]
    sp2 <- onDisk[[77]]
    expect_identical(sp1, sp2)

    ## by name.
    theN <- featureNames(inMem)[100]
    sp1 <- inMem[[theN]]
    sp2 <- onDisk[[theN]]
    expect_identical(sp1, sp2)

    sub1 <- inMem[1:20, ]
    sub2 <- onDisk[1:20, ]
    expect_true(all.equal(sub1, sub2))
    expect_true(all.equal(inMem[1], onDisk[1]))

    ## Test subsetting with file.
    sub1 <- filterAcquisitionNum(onDisk, n = 180:190)
})


############################################################
## validateOnDiskMSnExp
test_that("validateOnDiskMSnExp", {
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    onDisk <- readMSData2(f, verbose = FALSE)
    expect_true(validateOnDiskMSnExp(onDisk))
    ## Now modify the fData slightly.
    fd <- featureData(onDisk)
    fd$lowMZ[13] <- fd$lowMZ[13] + 3
    onDisk@featureData <- fd
    expect_error(validateOnDiskMSnExp(onDisk))
    suppressWarnings(
        expect_true(validateOnDiskMSnExp(filterFile(odmse, 1)))
    )
    ## That will cause an error; eventually there has been some data
    ## manipulations in that file?
    ## validateOnDiskMSnExp(filterFile(odmse, 2))
})


############################################################
## - other stuff -
## Compare the performacen of the C-contructor against the "standard" R constructor.
.compareCconstructorPerformance <- function(){
    featDat <- fData(odmse)
    featDat <- featDat[featDat$fileIdx == 1, ]
    ## Get all spectra from one file using the C-constructor
    system.time(
        spC <- MSnbase:::.applyFun2SpectraOfFile(featDat, filenames=fileNames(odmse))
    ) ## 3.7 sec; 4.2 sec
    ## Get all spectra from one file using the "new" constructor
    system.time(
        spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmse))
    ) ## 19 sec.; 24.2 sec
    system.time(
        spM <- MSnbase:::.applyFun2SpectraOfFileMulti(featDat, filenames=fileNames(odmse))
    ) ## 19 sec.; 3.1 sec
    expect_identical(spC, spR)
    expect_identical(spC, spM)
    ## Construct all of the spectra in one go...

    featDat <- fData(odmse)
    featDat <- featDat[featDat$fileIdx == 1, ]

    fileh <- mzR::openMSfile(fileNames(odmse)[1])
    ## Reading all of the data in "one go".
    allSpect <- mzR::peaks(fileh, featDat$spIdx)
    mzR::close(fileh)

    nValues <- lengths(allSpect) / 2
    allSpect <- do.call(rbind, allSpect)
    Test <- MSnbase:::Spectra1(peaksCount=featDat$peaksCount,
                               rt=featDat$retentionTime,
                               acquisitionNum=featDat$acquisitionNum,
                               tic=featDat$totIonCurrent,
                               mz=allSpect[, 1], intensity=allSpect[, 2],
                               centroided=featDat$centroided,
                               fromFile=rep(1, length(nValues)),
                               nvalues=nValues)
    names(Test) <- rownames(featDat)
    ## Have to change some stuff:
    ## o scanIndex is numeric() for the "standard" constructor
    ## o polarity is numeric().
    Test <- lapply(Test, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    ## This should be identical to the spectra from mse
    mseSpec <- spectra(mse)
    mseSpec <- mseSpec[fromFile(mse) == 1]
    expect_identical(mseSpec, Test)
}
