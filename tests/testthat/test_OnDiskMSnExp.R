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
## All the same with removePeaks.
mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram",
                          removePeaks=10000, clean=TRUE)
odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk",
                            removePeaks=10000, clean=TRUE)



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

test_that("rtime for OnDiskMSnExp", {
    rt <- rtime(mse)
    rt2 <- rtime(odmse)
    expect_identical(rt, rt2)
})

test_that("polarity for OnDiskMSnExp", {
    ##pol <- polarity(mse)
    pol2 <- polarity(odmse)
    ##expect_identical(pol, pol2)
})

test_that("tic for OnDiskMSnExp", {
    tc <- tic(mse)
    tc2 <- tic(odmse)
    expect_identical(tc, tc2)
})

test_that("ionCount for OnDiskMSnExp", {
    ic <- ionCount(mse)
    ic2 <- ionCount(odmse)
    expect_identical(ic, ic2)
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
    mseRemPeaks2 <- readMSData(files=mzf, msLevel=1, backend="ram",
                               removePeaks=100000)
    odmseRemPeaks2 <- readMSData(files=mzf, msLevel=1, backend="disk",
                                 removePeaks=100000)
    system.time(
        pk <- peaksCount(mseRemPeaks2)
    )  ## 0.046
    system.time(
        pk2 <- peaksCount(odmseRemPeaks2)
    )  ## 0 secs.
    expect_identical(pk, pk2)

    system.time(
        pk <- peaksCount(mseRemPeaks)
    )  ## 0.046
    ## Apply on spectra.
    system.time(
        pk2 <- peaksCount(odmseRemPeaks)
    )  ## 19 secs.
    expect_identical(pk, pk2)

    ## Peaks count with scans defined.
    system.time(
        ps3 <- peaksCount(mseRemPeaks, scans=c(1, 5, 9))
    ) ## 0.006
    system.time(
        ps4 <- peaksCount(odmseRemPeaks, scans=c(1, 5, 9))
    ) ## 0.05
    expect_identical(ps3, ps4)

    ## The same than pk2?
    expect_identical(ps4, pk2[c(1, 5, 9)])
})

test_that("compare spectra call", {
    system.time(
        spct <- spectra(mse)
    )  ## 0.005 sec.
    system.time(
        spct1 <- spectra(odmse)
    )  ## 7.1 sec.
    ## Polarity and scanIndex will not match.
    spct1 <- lapply(spct1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spct, spct1)

    system.time(
        spct <- spectra(mseRemPeaks)
    )  ## 0.005 sec.
    system.time(
        spct1 <- spectra(odmseRemPeaks)
    )  ## 18.8 sec.
    ## Polarity and scanIndex will not match.
    spct1 <- lapply(spct1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spct, spct1)

    ## Spectra for subsets.
    ## Note that ordering is NOT considered.
    subs <- c(4, 6, 13, 45, 3)
    system.time(
        spSub1 <- spectra(odmse, scans=subs)
    )  ## 0.068 sec
    ## Check if we really have the expected scans.
    ## Compare to manually extracted ones.
    spSub <- spectra(mse)[subs]
    spSub1 <- lapply(spSub1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spSub, spSub1)
    ## With the processing steps...
    system.time(
        spSub1 <- spectra(odmseRemPeaks, scans=subs)
    )  ## 0.088 sec
    spSub <- spectra(mseRemPeaks)[subs]
    spSub1 <- lapply(spSub1, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spSub, spSub1)

    ## Extract a single spectrum.
    sp1 <- spectra(mse)[[1]]
    sp2 <- spectra(odmse, scans=1)[[1]]
    sp2@polarity <- integer()
    sp2@scanIndex <- integer()
    expect_identical(sp1, sp2)
})

test_that("assayData on an OnDiskMSnExp", {
    env1 <- assayData(mse)
    env2 <- assayData(odmse)
    ## Can not compare directly, because polarity and scanIndex differ.
    expect_identical(ls(env1), ls(env2))
    env1 <- assayData(mseRemPeaks)
    env2 <- assayData(odmseRemPeaks)
    a <- unlist(eapply(env1, peaksCount))
    b <- unlist(eapply(env2, peaksCount))
    expect_identical(a, b[names(a)])
})

## removePeaks
test_that("removePeaks method for OnDiskMSnExp", {
    odmse2 <- removePeaks(odmse, t=10000)
    mse2 <- removePeaks(mse, t=10000)
    spOdmse2 <- spectra(odmse2)
    spMse2 <- spectra(mse2)
    spOdmse2 <- lapply(spOdmse2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spMse2, spOdmse2)
})

## clean
test_that("clean method for OnDiskMSnExp", {
    odmse2 <- removePeaks(odmse)
    mse2 <- removePeaks(mse)
    odmse2 <- clean(odmse2)
    mse2 <- clean(mse2)

    spOdmse2 <- spectra(odmse2)
    spMse2 <- spectra(mse2)
    spOdmse2 <- lapply(spOdmse2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(spMse2, spOdmse2)
    ##expect_identical(spectra(mse), spOdmse2)
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





