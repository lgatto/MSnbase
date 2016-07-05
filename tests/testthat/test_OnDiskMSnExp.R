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
mse <- readMSData(files = mzf, msLevel = 1, centroided = TRUE, verbose = FALSE)
## Load the data as OnDiskMSnExp.
odmse <- readMSData2(files = mzf, msLevel = 1, centroided=TRUE, verbose = FALSE)

## All the same with removePeaks.
mseRemPeaks <- readMSData(files = mzf, msLevel = 1, removePeaks = 10000,
                          clean = TRUE, verbose = FALSE)
odmseRemPeaks <- readMSData2(files = mzf, msLevel = 1, removePeaks = 10000,
                             clean = TRUE, verbose = FALSE)

############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
    expect_identical(as.character(class(mse)), "MSnExp")
    expect_identical(as.character(class(odmse)), "OnDiskMSnExp")
})

############################################################
## compare MSnExp against OnDiskMSnExp
test_that("Compare MS1 MSnExp and OnDiskMSnExp content", {
    ## Compare spectra values.
    expect_true(all.equal(mse, odmse))
    ## fromFile
    expect_identical(fromFile(mse), fromFile(odmse))
    ## msLevel
    expect_identical(msLevel(mse), msLevel(odmse))
})

############################################################
## acquisitionNum
test_that("compare acquisitionNum", {
    ## acquisitionNum
    expect_identical(acquisitionNum(mse), acquisitionNum(odmse))
})

############################################################
## scanIndex
test_that("compare scanIndex", {
    ## scanIndex; point is that scanIndex on an MSnExp will return 0, as it is not
    ## set in the spectra (same as acquisitionNum?)
    expect_identical(acquisitionNum(mse), scanIndex(odmse))
})

############################################################
## centroided
test_that("compare centroided", {
    ## centroided.
    expect_identical(centroided(mse), centroided(odmse))
    ## Setting stuff
    centroided(mse) <- FALSE  ## Takes quite some time.
    centroided(odmse) <- FALSE
    expect_identical(centroided(mse), centroided(odmse))
    ## Check error
    expect_that(centroided(odmse) <- c(TRUE, FALSE, TRUE), throws_error())
})

############################################################
## rtime
test_that("rtime for OnDiskMSnExp", {
    rt <- rtime(mse)
    rt2 <- rtime(odmse)
    expect_identical(rt, rt2)
})

############################################################
## polarity
test_that("polarity for OnDiskMSnExp", {
    pol <- polarity(mse)
    pol2 <- polarity(odmse)
    expect_identical(pol, pol2)
})

############################################################
## tic
test_that("ionCount for OnDiskMSnExp", {
    tc <- tic(mse)
    tc2 <- tic(odmse)
    expect_identical(tc, tc2)
})

############################################################
## ionCount
test_that("ionCount for OnDiskMSnExp", {
    system.time(
        ic <- ionCount(mse)
    ) ## 0.058
    system.time(
        ic2 <- ionCount(odmse)
    ) ## 5 sec.
    expect_identical(ic, ic2)
})

############################################################
## peaksCount
test_that("compare peaksCount", {
    ## Trivial case, without any processing steps.
    ## peaksCount
    system.time(
        pk <- peaksCount(mse)
    )  ## 0.049
    system.time(
        pk2 <- peaksCount(odmse)
    )  ## 0.002
    expect_identical(pk, pk2)
})

############################################################
## Compare cleaned data.
## o spectra.
## o ionCount.
## o tic
## o peaksCount.
test_that("Compare cleaned MSnExp and OnDiskMSnExp", {
    mseCleaned <- clean(mse)
    odmseCleaned <- clean(odmse)
    ## o spectra:
    expect_true(all.equal(mseCleaned, odmseCleaned))

    ## o ionCount.
    expect_identical(ionCount(mseCleaned), ionCount(odmseCleaned))

    ## o tic.
    expect_identical(tic(mseCleaned), tic(odmseCleaned))

    ## o peaksCount
    expect_identical(peaksCount(mseCleaned), peaksCount(odmseCleaned))

})

############################################################
## Compare removePeaks data
## o spectra.
## o ionCount.
## o tic.
## o peaksCount.
test_that("Compare removePeaks results between MSnExp and OnDiskMSnExp", {
    mseRemP <- removePeaks(mse, t=1000)
    odmseRemP <- removePeaks(odmse, t=1000)
    ## o spectra:
    expect_true(all.equal(mseRemP, odmseRemP))

    ## o ionCount.
    expect_identical(ionCount(mseRemP), ionCount(odmseRemP))

    ## o tic.
    expect_identical(tic(mseRemP), tic(odmseRemP))

    ## o peaksCount
    expect_identical(peaksCount(mseRemP), peaksCount(odmseRemP))
})

############################################################
## Compare removePeaks and cleaned data
## o spectra.
## o ionCount.
## o tic.
## o peaksCount.
test_that("Compare removePeaks and clean results between MSnExp and OnDiskMSnExp", {
    ## o spectra:
    expect_true(all.equal(mseRemPeaks, odmseRemPeaks))

    ## o ionCount.
    expect_identical(ionCount(mseRemPeaks), ionCount(odmseRemPeaks))

    ## o tic.
    expect_identical(tic(mseRemPeaks), tic(odmseRemPeaks))

    ## o peaksCount
    expect_identical(peaksCount(mseRemPeaks), peaksCount(odmseRemPeaks))
})

############################################################
## assayData
test_that("assayData on an OnDiskMSnExp", {
    env1 <- assayData(mse)
    env2 <- assayData(odmse)
    expect_equal(env1, env2)

    env1 <- assayData(mseRemPeaks)
    env2 <- assayData(odmseRemPeaks)
    expect_equal(env1, env2)
})

############################################################
## intensity
test_that("intensity on an OnDiskMSnExp", {
    ints <- intensity(mse)
    system.time(
        ints2 <- intensity(odmse)
    ) ## 5.5 sec
    expect_identical(ints, ints2)

    ints <- intensity(mseRemPeaks)
    system.time(
        ints2 <- intensity(odmseRemPeaks)
    ) ## 18.3 sec
    expect_identical(ints, ints2)
})

############################################################
## mz
test_that("mz on an OnDiskMSnExp", {
    mzs <- mz(mse)
    system.time(
        mzs2 <- mz(odmse)
    ) ## 5.5 sec
    expect_identical(mzs, mzs2)

    mzs <- mz(mseRemPeaks)
    system.time(
        mzs2 <- mz(odmseRemPeaks)
    ) ## 18 sec
    expect_identical(mzs, mzs2)
})

############################################################
## [[
test_that("[[ for OnDiskMSnExp", {
    sp1 <- mse[[77]]
    sp2 <- odmse[[77]]
    expect_identical(sp1, sp2)

    ## by name.
    theN <- featureNames(mse)[100]
    sp1 <- mse[[theN]]
    sp2 <- odmse[[theN]]
    expect_identical(sp1, sp2)
})

############################################################
## [
test_that("[ for OnDiskMSnExp", {
    ## subset by row (i)
    sub1 <- mse[1:20, ]
    sub2 <- odmse[1:20, ]
    expect_true(all.equal(sub1, sub2))
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
