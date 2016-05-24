context("OnDiskMSnExp class")

############################################################
## Load the required data files.
.getMzMLFiles <- function(){
    ## Return the mzML files, the ones from the XXX package, or if run
    ## locally, some of my test files.
    HOST <- unlist(strsplit(system("hostname", intern=TRUE), split=".",
                            perl=FALSE, fixed=TRUE))[1]
    if(HOST == "macbookjo"){
        mzfiles <- dir("/Users/jo/R-workspaces/EURAC/2016/2016-04-21-PolarMetabolom/data/mzML/",
                       pattern="POS_C_O", full.names=TRUE)
    }else{
        require(msdata)
        mzfiles <- c(system.file("microtofq/MM14.mzML", package="msdata"),
                     system.file("microtofq/MM8.mzML", package="msdata"))
    }
    return(mzfiles)
}
mzf <- .getMzMLFiles()[1:2]
## Load the data as an MSnExp into memory.
mse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="ram")
## Load the data as OnDiskMSnExp.
odmse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="disk")



############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
    expect_identical(as.character(class(mse)), "MSnExp")
    expect_identical(as.character(class(odmse)), "OnDiskMSnExp")
})

test_that("compare basic contents", {
    ## Check if we get the same data! MSnExp will retrieve that from the spectra,
    ## OnDiskMSnExp from featureData.
    ## fromFile
    expect_identical(fromFile(mse), fromFile(odmse))

    ## msLevel
    expect_identical(msLevel(mse), msLevel(odmse))

    ## header
    hd <- header(mse)
    odhd <- header(odmse)
    commonCols <- intersect(colnames(hd), colnames(odhd))
    expect_identical(hd[, commonCols], odhd[, commonCols])

    ## length
    expect_identical(length(mse), length(odmse))

})

test_that("compare acquisitionNum", {
    ## acquisitionNum
    expect_identical(acquisitionNum(mse), acquisitionNum(odmse))
})

test_that("compare scanIndex", {
    ## scanIndex; point is that scanIndex on an MSnExp will return 0, as it is not
    ## set in the spectra (same as acquisitionNum?)
    expect_identical(acquisitionNum(mse), scanIndex(odmse))
})

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

    ## Have a processing step "removePeaks"
    mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram",
                              removePeaks=100000)
    odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk",
                                removePeaks=100000)
    system.time(
        pk <- peaksCount(mseRemPeaks)
    )  ## 0.046
    system.time(
        pk2 <- peaksCount(odmseRemPeaks)
    )  ## 0 secs.
    expect_identical(pk, pk2)

    ## Now including also "clean", i.e. this means we'll have to recalculate the
    ## peaks count on-the-fly.
    mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram",
                              removePeaks=100000, clean=TRUE)
    odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk",
                                removePeaks=100000, clean=TRUE)
    system.time(
        pk <- peaksCount(mseRemPeaks)
    )  ## 0.046
    ## Apply on spectra.
    system.time(
        pk2 <- peaksCount(odmseRemPeaks)
    )  ## 26 secs.
    expect_identical(pk, pk2)
    system.time(
        pk3 <- peaksCount(odmseRemPeaks, method=2)
    )  ## 25 secs.
    expect_identical(pk, pk3)
})

test_that("compare spectra call", {
    system.time(
        spct <- spectra(mse)
    )  ## 0.005 sec.
    system.time(
        spct1 <- spectra(odmse)
    )  ## 7.1 sec.
    system.time(
        spct2 <- spectra(odmse, method=2)
    )  ## 5.4 sec.
    expect_identical(spct1, spct2)
    ## Polarity and scanIndex will not match.
    spct1 <- lapply(spct1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    spct2 <- lapply(spct2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spct, spct1)
    expect_identical(spct, spct2)

    ## All the same with removePeaks.
    mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram",
                              removePeaks=10000, clean=TRUE)
    odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk",
                                removePeaks=10000, clean=TRUE)
    system.time(
        spct <- spectra(mseRemPeaks)
    )  ## 0.005 sec.
    system.time(
        spct1 <- spectra(odmseRemPeaks)
    )  ## 28.8 sec.
    system.time(
        spct2 <- spectra(odmseRemPeaks, method=2)
    )  ## 25.6 sec.
    ## Polarity and scanIndex will not match.
    expect_identical(spct1, spct2)
    spct1 <- lapply(spct1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    spct2 <- lapply(spct2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spct, spct1)
    expect_identical(spct, spct2)

    ## Spectra for subsets.
    ## Note that ordering is NOT considered.
    subs <- c(4, 6, 13, 45, 3)
    system.time(
        spSub1 <- spectra(odmse, scans=subs, method=1)
    )  ## 0.068 sec
    system.time(
        spSub2 <- spectra(odmse, scans=subs, method=2)
    )  ## 0.052
    expect_identical(spSub1, spSub2)
    ## Check if we really have the expected scans.
    ## Compare to manually extracted ones.
    spSub <- spectra(mse)[acquisitionNum(mse) %in% subs]
    spSub1 <- lapply(spSub1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    spSub2 <- lapply(spSub2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spSub, spSub1)
    expect_identical(spSub, spSub2)
    expect_identical(spSub2, spSub1)
    ## With the processing steps...
    system.time(
        spSub1 <- spectra(odmseRemPeaks, scans=subs, method=1)
    )  ## 0.088 sec
    system.time(
        spSub2 <- spectra(odmseRemPeaks, scans=subs, method=2)
    )  ## 0.082
    expect_identical(spSub1, spSub2)
    spSub <- spectra(mseRemPeaks)[acquisitionNum(mseRemPeaks) %in% subs]
    spSub1 <- lapply(spSub1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    spSub2 <- lapply(spSub2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spSub, spSub1)
    expect_identical(spSub, spSub2)
    expect_identical(spSub2, spSub1)
})


## Compare the performacen of the C-contructor against the "standard" R constructor.
.compareCconstructorPerformance <- function(){
    featDat <- fData(odmse)
    featDat <- featDat[featDat$fileIdx == 1, ]
    ## Get all spectra from one file using the C-constructor
    system.time(
        spC <- MSnbase:::.applyFun2SpectraOfFile(featDat, filenames=fileNames(odmse))
    ) ## 3.7 sec.
    ## Get all spectra from one file using the "new" constructor
    system.time(
        spR <- MSnbase:::.applyFun2SpectraOfFileSlow(featDat, filenames=fileNames(odmse))
    ) ## 19 sec.
    expect_identical(spC, spR)
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





