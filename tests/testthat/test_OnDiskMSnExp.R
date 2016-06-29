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

############################################################
## compare MSnExp against OnDiskMSnExp
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
    expect_equal(hd[, commonCols], odhd[, commonCols])

    ## length
    expect_identical(length(mse), length(odmse))

})

############################################################
## header
test_that("header on OnDiskMSnExp", {
    system.time(
        hd1 <- header(mse)
    )
    system.time(
        hd2 <- header(odmse)
    )
    commonCols <- intersect(colnames(hd1), colnames(hd2))
    expect_equal(hd1[, commonCols], hd2[, commonCols])

    ## header with scans.
    system.time(
        hd1 <- header(mse, scans=1:300)
    ) ## 0.5
    system.time(
        hd2 <- header(odmse, scans=1:300)
    ) ## 0.3
    commonCols <- intersect(colnames(hd1), colnames(hd2))
    expect_equal(hd1[, commonCols], hd2[, commonCols])
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
    ##pol <- polarity(mse)
    pol2 <- polarity(odmse)
    ##expect_identical(pol, pol2)
})

############################################################
## tic
test_that("tic for OnDiskMSnExp", {
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

############################################################
## spectra
test_that("compare spectra call", {
    system.time(
        spct <- spectra(mse)
    )  ## 0.005 sec.
    system.time(
        spct1 <- spectra(odmse)
    )  ## 5.6 sec.
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
    )  ## 20 sec.
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
    )  ## 0.069 sec
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

    ## Extract by retention time limit.
    rtl <- c(4, 20)
    subs <- rtime(mse) >= rtl[1] & rtime(mse) <= rtl[2]
    sp1 <- spectra(mse)[subs]

    sp2 <- spectra(odmse, rtlim=rtl)
    sp2 <- lapply(sp2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(sp1, sp2)
})

############################################################
## assayData
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

############################################################
## intensity
test_that("intensity on an OnDiskMSnExp", {
    ints <- intensity(mse)
    system.time(
        ints2 <- intensity(odmse)
    ) ## 5.5 sec
    expect_identical(ints, ints2)

    ## Now what if we've got some processings...
    mse2 <- removePeaks(mse)
    odmse2 <- removePeaks(odmse)
    ints <- intensity(mse2)
    system.time(
        ints2 <- intensity(odmse2)
    ) ## 6.9 sec
    expect_identical(ints, ints2)

    ## Now what if we've got some processings...
    mse2 <- clean(mse2)
    odmse2 <- clean(odmse2)
    ints <- intensity(mse2)
    system.time(
        ints2 <- intensity(odmse2)
    ) ## 49.9 sec; eventually we might speed up the clean at some point!
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

    ## Now what if we've got some processings...
    mse2 <- removePeaks(mse)
    odmse2 <- removePeaks(odmse)
    mzs <- mz(mse2)
    system.time(
        mzs2 <- mz(odmse2)
    ) ## 6.9 sec
    expect_identical(mzs, mzs2)

    ## Now what if we've got some processings...
    mse2 <- clean(mse2)
    odmse2 <- clean(odmse2)
    mzs <- mz(mse2)
    system.time(
        mz2 <- mz(odmse2)
    ) ## 49.9 sec; eventually we might speed up the clean at some point!
    expect_identical(mzs, mz2)
})

############################################################
## [[
test_that("[[ for OnDiskMSnExp", {
    sp1 <- mse[[77]]
    sp2 <- odmse[[77]]
    sp2@polarity <- integer()
    sp2@scanIndex <- integer()
    expect_identical(sp1, sp2)

    ## by name.
    theN <- featureNames(mse)[100]
    sp1 <- mse[[theN]]
    sp2 <- odmse[[theN]]
    sp2@polarity <- integer()
    sp2@scanIndex <- integer()
    expect_identical(sp1, sp2)
})

############################################################
## [
test_that("[ for OnDiskMSnExp", {
    ## subset by row (i)
    sub1 <- mse[1:20, ]
    sub2 <- odmse[1:20, ]
    expect_identical(featureNames(sub1), featureNames(sub2))
    sp1 <- spectra(sub1)
    sp2 <- spectra(sub2)
    sp2 <- lapply(sp2, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
        return(z)
    })
    expect_identical(sp1, sp2)

    ## subset by sample (j)
    sub2 <- odmse[, 2]
    ## featureData I expect to have only fileIndex 1
    expect_equal(unique(fData(sub2)$fileIdx), 1)
    expect_identical(fileNames(sub2), fileNames(odmse)[2])
    expect_identical(unique(fileNames(sub2)[fData(sub2)$fileIdx]), fileNames(odmse)[2])
    ## compare the featureData
    fidxCol <- which(colnames(fData(sub2)) == "fileIdx")
    expect_identical(fData(sub2)[, -fidxCol], fData(odmse)[fData(odmse)$fileIdx == 2, -fidxCol])
    ## compare spectra
    sp1 <- spectra(sub2)
    sp2 <- spectra(odmse, scans=fData(odmse)$fileIdx == 2)
    ## Fix fromFile.
    sp2 <- lapply(sp2, function(z){
        z@fromFile <- 1L
        return(z)
    })
    expect_identical(sp1, sp2)

    ## subset by i and j
    sub1 <- odmse[1:20, 2]
    sub2 <- odmse[1:20, ]
    sp1 <- spectra(sub1)
    sp2 <- spectra(sub2, scans=fData(sub2)$fileIdx == 2)
    ## Fix fromFile.
    sp2 <- lapply(sp2, function(z){
        z@fromFile <- 1L
        return(z)
    })
    expect_identical(sp1, sp2)
})

############################################################
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

############################################################
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





