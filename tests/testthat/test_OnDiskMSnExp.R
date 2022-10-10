context("OnDiskMSnExp class")

inMem <- microtofq_in_mem_ms1
onDisk <- microtofq_on_disk_ms1
multiMsInMem1 <- tmt_im_ms1_sub
multiMsInMem2 <- tmt_im_ms2_sub
multiMsOnDisk <- tmt_od_sub
centroided(inMem) <- TRUE
centroided(onDisk) <- TRUE
centroided(multiMsInMem1) <- TRUE
centroided(multiMsInMem2) <- TRUE
centroided(multiMsOnDisk) <- TRUE

############################################################
## validateOnDiskMSnExp
test_that("validateOnDiskMSnExp", {
    onDisk2 <- multiMsOnDisk
    expect_true(validateOnDiskMSnExp(onDisk2))
    ## Now modify the fData slightly.
    fd <- featureData(onDisk2)
    fd$lowMZ[13] <- fd$lowMZ[13] + 3
    onDisk2@featureData <- fd
    expect_error(validateOnDiskMSnExp(onDisk2))
    expect_true(validateOnDiskMSnExp(filterMsLevel(onDisk2, 2)))
})

test_that("Empty validateOnDiskMSnExp", {
    x <- filterMsLevel(onDisk, 2) ## empty
    expect_identical(length(x), 0L)
    expect_true(validObject(x))
    expect_true(validateOnDiskMSnExp(x))
})

test_that("Warning validateOnDiskMSnExp", {
    expect_warning(val <- validateOnDiskMSnExp(onDisk))
    expect_true(val)
})

############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
    expect_identical(as.character(class(inMem)), "MSnExp")
    expect_identical(as.character(class(onDisk)), "OnDiskMSnExp")
    expect_true(validObject(onDisk))
})

test_that("Coercion to MSnExp", {
    x <- tmt_erwinia_on_disk
    y <- tmt_erwinia_in_mem_ms2
    expect_error(as(x, "MSnExp"))
    x <- filterMsLevel(x, msLevel = 2)
    expect_true(all.equal(x, y))
    ## feature names are different
    featureNames(x) <- featureNames(y)
    x2 <- as(x, "MSnExp")
    ## expected to be different: processingData, fData, .cache
    expect_identical(spectra(x2), spectra(y))
    expect_identical(experimentData(x2), experimentData(y))
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
    ## o tic: Note: have to calculate as one of the two doesn't provide the
    ##        initial values.
    expect_identical(tic(inMem), tic(onDisk, initial = FALSE))
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
    expect_identical(tic(inMemCleaned), tic(onDiskCleaned, initial = FALSE))
    expect_identical(peaksCount(inMemCleaned), peaksCount(onDiskCleaned))
    ## o removePeaks
    inMemRemPeaks <- removePeaks(inMem, t = 1000)
    onDiskRemPeaks <- removePeaks(onDisk, t = 1000)
    expect_true(all.equal(inMemRemPeaks, onDiskRemPeaks))
    expect_identical(ionCount(inMemRemPeaks), ionCount(onDiskRemPeaks))
    expect_identical(tic(inMemRemPeaks), tic(onDiskRemPeaks, initial = FALSE))
    expect_identical(peaksCount(inMemRemPeaks), peaksCount(onDiskRemPeaks))
    ## o removePeaks and clean
    inMemRemPeaksCleaned <- clean(inMemRemPeaks)
    onDiskRemPeaksCleaned <- clean(onDiskRemPeaks)
    expect_true(all.equal(inMemRemPeaksCleaned, onDiskRemPeaksCleaned))
    expect_identical(ionCount(inMemRemPeaksCleaned),
                     ionCount(onDiskRemPeaksCleaned))
    expect_identical(tic(inMemRemPeaksCleaned), tic(onDiskRemPeaksCleaned,
                                                    initial = FALSE))
    expect_identical(peaksCount(inMemRemPeaksCleaned),
                     peaksCount(onDiskRemPeaksCleaned))
    ## compare assayData, intensity and mz,
    expect_equal(assayData(inMemRemPeaksCleaned),
                 assayData(onDiskRemPeaksCleaned))
    expect_equal(intensity(inMemRemPeaksCleaned),
                 intensity(onDiskRemPeaksCleaned))
    expect_equal(mz(inMemRemPeaksCleaned), mz(onDiskRemPeaksCleaned))
})

test_that("clean on OnDiskMSnExp with different MS levels", {
    ## o Tests on MSnExp
    multiMsInMem1_cleaned <- clean(multiMsInMem1)
    expect_true(sum(unlist(intensity(multiMsInMem1_cleaned)) == 0) <
                sum(unlist(intensity(multiMsInMem1)) == 0))
    ## o Tests on OnDiskMSnExp and comparison with MSnExp.
    multiMsOnDisk_cleaned <- clean(multiMsOnDisk)
    expect_true(sum(unlist(intensity(multiMsOnDisk_cleaned)) == 0) <
                sum(unlist(intensity(multiMsOnDisk)) == 0))
    ##   Compare with MSnExp
    expect_true(all.equal(multiMsInMem1_cleaned,
                          filterMsLevel(multiMsOnDisk_cleaned, msLevel. = 1)))
    ##   Just cleaning MS 1.
    multiMsOnDisk_cleaned_1 <- clean(multiMsOnDisk, msLevel. = 1)
    expect_true(all.equal(multiMsOnDisk_cleaned, multiMsOnDisk_cleaned_1))
    ##   Just cleaning MS 2; won't do much at all.
    multiMsOnDisk_cleaned_2 <- clean(multiMsOnDisk, msLevel. = 2)
    expect_true(all.equal(multiMsOnDisk, multiMsOnDisk_cleaned_2))
    ##   Same with msLevel. 4
    multiMsOnDisk_cleaned_4 <- clean(multiMsOnDisk, msLevel. = 4)
    expect_true(all.equal(multiMsOnDisk, multiMsOnDisk_cleaned_4))
})

test_that("removePeaks on OnDiskMSnExp with different MS levels", {
    ## o Tests on MSnExp
    multiMsInMem1_rem <- removePeaks(multiMsInMem1)
    expect_true(sum(unlist(intensity(multiMsInMem1_rem)) == 0) >
                sum(unlist(intensity(multiMsInMem1)) == 0))
    multiMsInMem2_rem <- removePeaks(multiMsInMem2)
    expect_true(sum(unlist(intensity(multiMsInMem2_rem)) == 0) >
                sum(unlist(intensity(multiMsInMem2)) == 0))
    ## o Tests on OnDiskMSnExp and comparison with MSnExp.
    multiMsOnDisk_rem <- removePeaks(multiMsOnDisk)
    expect_true(sum(unlist(intensity(multiMsOnDisk_rem)) == 0) >
                sum(unlist(intensity(multiMsOnDisk)) == 0))
    ##   Compare with MSnExp
    expect_true(all.equal(multiMsInMem1_rem,
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 1)))
    expect_true(all.equal(multiMsInMem2_rem,
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 2)))
    ##   Just processing MS 1.
    multiMsOnDisk_rem_1 <- removePeaks(multiMsOnDisk, msLevel. = 1)
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 1),
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 1)))
    spects1 <- spectra(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 2))
    spects2 <- spectra(filterMsLevel(multiMsOnDisk, msLevel. = 2))
    expect_identical(spects1, spects2)
    ## expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 2),
    ##                       filterMsLevel(multiMsOnDisk, msLevel. = 2)))
    ##   Just processing MS 2.
    multiMsOnDisk_rem_2 <- removePeaks(multiMsOnDisk, msLevel. = 2)
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_2, msLevel. = 2),
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 2)))
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_2, msLevel. = 1),
                          filterMsLevel(multiMsOnDisk, msLevel. = 1)))
})


############################################################
## bin
test_that("bin on OnDiskMSnExp", {
    ## o On a single-file multi MS-level file.
    multiMsInMem1_bin <- bin(multiMsInMem1, verbose = FALSE)
    multiMsInMem2_bin <- bin(multiMsInMem2, verbose = FALSE)
    ##   bin on MS1 level only
    multiMsOnDisk_bin_1 <- bin(multiMsOnDisk, msLevel. = 1)
    ##   Results should be the same.
    expect_true(all.equal(multiMsInMem1_bin,
                          filterMsLevel(multiMsOnDisk_bin_1, msLevel. = 1)))
    ##   bin on all levels.
    multiMsOnDisk_bin <- bin(multiMsOnDisk)
    ##   Results can not be the same, since the mz range was different for
    ##   the bin.
    expect_true(is(all.equal(
        multiMsInMem1_bin,
        filterMsLevel(multiMsOnDisk_bin, msLevel. = 1)
    ), "character"))
    ##  bin on MS2 level only
    multiMsOnDisk_bin_2 <- bin(multiMsOnDisk, msLevel. = 2)
    ##   Results should be the same.
    expect_true(all.equal(multiMsInMem2_bin,
                          filterMsLevel(multiMsOnDisk_bin_2, msLevel. = 2)))

    ## o On multiple files.
    inMem_bin <- bin(inMem, verbose = FALSE)
    onDisk_bin <- bin(onDisk)
    expect_true(all.equal(inMem_bin, onDisk_bin))
    ##   bin on MS 2 shouldn't do anything at all
    expect_warning(onDisk_bin <- bin(onDisk, msLevel. = 2))
    expect_true(all.equal(onDisk_bin, onDisk))
})


############################################################
## Test internal spectrapply method.
test_that("Test internal spectrapply function", {
    spl <- spectra(onDisk)
    ## Test Spectrum method:
    res1 <- lapply(spl, ionCount)
    res2 <- spectrapply(onDisk, ionCount)
    expect_identical(res1, res2)

    ## Plain function
    res1 <- lapply(spl, function(z) return(mean(mz(z))))
    res2 <- spectrapply(onDisk, function(z) {
        return(mean(mz(z)))
    })
    expect_identical(res1, res2)

    ## Additional arguments.
    res1 <- lapply(spl, function(z, int) {
        return(mean(mz(z)[intensity(z) > int]))
    }, int = 30)
    res2 <- spectrapply(onDisk, function(z, int) {
        return(mean(mz(z)[intensity(z) > int]))
    }, int = 30)
    expect_identical(res1, res2)
})

############################################################
## Test that the new sorting by acquisitionNum and extraction
## by spIdx works (described in issue #118)
## We're comparing spectra extracted by that between an
## OnDiskMSnExp and an MSnExp.
test_that("Test sorting by acquisitionNum", {
    sp1 <- inMem[[13]]
    sp2 <- onDisk[[13]]
    expect_identical(sp1, sp2)
    sp1 <- inMem[[22]]
    sp2 <- onDisk[[22]]
    expect_identical(sp1, sp2)
    ## Same using multiMS
    onDisk1 <- filterMsLevel(multiMsOnDisk, msLevel. = 1L)
    sp1 <- multiMsInMem1[[7]]
    sp2 <- onDisk1[[7]]
    expect_identical(sp1, sp2)
    sp1 <- multiMsInMem1[[9]]
    sp2 <- onDisk1[[9]]
    expect_identical(sp1, sp2)
    onDisk2 <- filterMsLevel(multiMsOnDisk, msLevel. = 2L)
    sp1 <- multiMsInMem2[[13]]
    sp2 <- onDisk2[[13]]
    expect_identical(sp1, sp2)
    sp1 <- multiMsInMem2[[22]]
    sp2 <- onDisk2[[22]]
    expect_identical(sp1, sp2)
})

test_that("spectrapply,OnDiskMSnExp", {
    sps <- spectra(onDisk)
    sps_2 <- spectrapply(onDisk)
    expect_identical(sps, sps_2)
    ## apply a function.
    dfs <- spectrapply(onDisk, FUN = as, Class = "data.frame")
    dfs_2 <- lapply(sps, FUN = as, Class = "data.frame")
    expect_identical(dfs, dfs_2)
})

test_that("splitByFile,OnDiskMSnExp", {
    od <- microtofq_on_disk_ms1
    spl <- splitByFile(od, f = factor(c("a", "b")))
    expect_equal(pData(spl[[1]]), pData(filterFile(od, 1)))
    expect_equal(pData(spl[[2]]), pData(filterFile(od, 2)))
})

test_that("chromatogram,OnDiskMSnExp works", {
    library(msdata)
    mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
             system.file("microtofq/MM8.mzML", package = "msdata"))
    tmpd <- tempdir()
    file.copy(mzf[1], paste0(tmpd, "a.mzML"))
    file.copy(mzf[2], paste0(tmpd, "b.mzML"))
    mzf <- c(mzf, paste0(tmpd, c("a.mzML", "b.mzML")))

    onDisk <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE,
                         mode = "onDisk")

    ## Full rt range.
    mzr <- matrix(c(100, 120), nrow = 1)
    res <- MSnbase:::.extractMultipleChromatograms(onDisk, mz = mzr)
    flt <- filterMz(onDisk, mz = mzr[1, ])
    ints <- split(unlist(lapply(spectra(flt), function(z) sum(intensity(z)))),
                  fromFile(flt))
    expect_equal(ints[[1]], intensity(res[1, 1][[1]]))
    expect_equal(ints[[2]], intensity(res[1, 2][[1]]))
    expect_equal(split(rtime(flt), fromFile(flt))[[1]], rtime(res[1, 1][[1]]))
    expect_equal(split(rtime(flt), fromFile(flt))[[2]], rtime(res[1, 2][[1]]))

    ## Multiple mz ranges.
    mzr <- matrix(c(100, 120, 200, 220, 300, 320), nrow = 3, byrow = TRUE)
    rtr <- matrix(c(50, 300), nrow = 1)
    res <- MSnbase:::.extractMultipleChromatograms(onDisk, mz = mzr, rt = rtr)
    ## Check that the values for all ranges is within the specified ranges
    for (i in 1:nrow(mzr)) {
        expect_true(all(mz(res[i, 1][[1]]) >= mzr[i, 1] &
                        mz(res[i, 1][[1]]) <= mzr[i, 2]))
        expect_true(all(mz(res[i, 2][[1]]) >= mzr[i, 1] &
                        mz(res[i, 2][[1]]) <= mzr[i, 2]))
        expect_true(all(rtime(res[i, 1][[1]]) >= rtr[1, 1] &
                        rtime(res[i, 1][[1]]) <= rtr[1, 2]))
        expect_true(all(rtime(res[i, 2][[1]]) >= rtr[1, 1] &
                        rtime(res[i, 2][[1]]) <= rtr[1, 2]))
    }
    ## Check that values are correct.
    flt <- filterMz(filterRt(onDisk, rt = rtr[1, ]), mz = mzr[2, ])
    ints <- split(unlist(lapply(spectra(flt), function(z) sum(intensity(z)))),
                  fromFile(flt))
    expect_equal(ints[[1]], intensity(res[2, 1][[1]]))
    expect_equal(ints[[2]], intensity(res[2, 2][[1]]))
    expect_equal(split(rtime(flt), fromFile(flt))[[1]], rtime(res[2, 1][[1]]))
    expect_equal(split(rtime(flt), fromFile(flt))[[2]], rtime(res[2, 2][[1]]))

    ## Now with ranges for which we don't have values in one or the other.
    rtr <- matrix(c(280, 300, 20, 40), nrow = 2,
                  byrow = TRUE)  ## Only present in first, or 2nd file
    res <- chromatogram(onDisk, rt = rtr)
    expect_true(all(unlist(lapply(res, msLevel)) == 1))
    ## Check fromFile
    for (i in 1:ncol(res))
        expect_true(all(sapply(res[, i], fromFile) == i))
    expect_equal(length(res[2, 1]), 0)
    expect_equal(length(res[1, 2]), 0)
    ## Check rtime
    expect_true(all(rtime(res[1, 1]) >= rtr[1, 1] &
                    rtime(res[1, 1]) <= rtr[1, 2]))
    expect_true(all(rtime(res[2, 2]) >= rtr[2, 1] &
                    rtime(res[2, 2]) <= rtr[2, 2]))
    ## Check intensity
    flt <- filterRt(onDisk, rt = rtr[1, ])
    spctr <- split(spectra(flt), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[1, 1]))
    flt <- filterRt(onDisk, rt = rtr[2, ])
    spctr <- split(spectra(flt), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[2, 2]))

    ## Check chromatogram with non-present MS level
    expect_warning(tmp <- chromatogram(onDisk, rt = rtr, msLevel = 2L))
    expect_equal(nrow(tmp), 0)
    tmp <- chromatogram(onDisk, rt = rtr, msLevel = 1:10)
    expect_equal(tmp, res)

    res <- MSnbase:::.extractMultipleChromatograms(onDisk, rt = rtr,
                                                   msLevel = 1:5)
    colnames(res) <- basename(fileNames(onDisk))
    res <- as(res, "MChromatograms")
    expect_true(validObject(res))
    ## pData(res) <- pData(onDisk)
    ## fData(res) <- fData(tmp)
    ## expect_equal(tmp, res)
})

## Test the two versions that could/might be called by the
## spectrapply,OnDiskMSnExp method. Each has its own pros and cons and cases
## in which it outperforms the other function.
test_that("low memory spectrapply function works", {
    fl <- system.file("lockmass/LockMass_test.mzXML", package = "msdata")
    fh <- mzR::openMSfile(fl)
    hdr <- mzR::header(fh)
    mzR::close(fh)
    fData <- hdr
    fData$spIdx <- hdr$seqNum
    fData$fileIdx <- 1L
    fData$smoothed <- FALSE
    fData$centroided <- TRUE

    fastLoad <- FALSE

    expect_equal(
        MSnbase:::.applyFun2SpectraOfFileMulti(fData, filenames = fl,
                                               fastLoad = fastLoad),
        MSnbase:::.applyFun2IndividualSpectraOfFile(fData,
                                                    filenames = fl,
                                                    fastLoad = fastLoad))
    fd <- fData[c(4, 8, 32, 123), ]
    expect_equal(
        MSnbase:::.applyFun2SpectraOfFileMulti(fd, filenames = fl,
                                               fastLoad = fastLoad),
        MSnbase:::.applyFun2IndividualSpectraOfFile(fd,
                                                    filenames = fl,
                                                    fastLoad = fastLoad))
    ## With an function to apply.
    expect_equal(
        MSnbase:::.applyFun2SpectraOfFileMulti(fd, filenames = fl,
                                               fastLoad = fastLoad,
                                               APPLYFUN = mz),
        MSnbase:::.applyFun2IndividualSpectraOfFile(fd,
                                                    filenames = fl,
                                                    fastLoad = fastLoad,
                                                    APPLYFUN = mz))
})

test_that("isolationWindowLowerMz,isolationWindowUpperMz,OnDiskMSnExp works", {
    mz_low <- isolationWindowLowerMz(tmt_od_ms2_sub)
    mz_high <- isolationWindowUpperMz(tmt_od_ms2_sub)
    expect_true(all(mz_low < mz_high))
    expect_true(all(precursorMz(tmt_od_ms2_sub) >= mz_low))
    expect_true(all(precursorMz(tmt_od_ms2_sub) <= mz_high))

    mz_low <- isolationWindowLowerMz(sciex)
    mz_high <- isolationWindowUpperMz(sciex)
    expect_true(all(is.na(mz_low)))
    expect_true(all(is.na(mz_high)))
})

test_that("combineSpectra,MSnExp works with OnDiskMSnExp", {
    res <- combineSpectra(filterRt(sciex, c(10, 20)))
    expect_true(is(res, "MSnExp"))
    expect_true(length(res) == 2)
})

test_that(".on_disk_split_by_file works", {
    res <- .on_disk_split_by_file(sciex)
    expect_equal(length(res), 2)
    expect_equal(featureData(res[[1]]), featureData(filterFile(sciex, 1L)))
    expect_equal(featureData(res[[2]]), featureData(filterFile(sciex, 2L)))
    expect_equal(phenoData(res[[1]]), phenoData(filterFile(sciex, 1L)))
    expect_equal(phenoData(res[[2]]), phenoData(filterFile(sciex, 2L)))
})
