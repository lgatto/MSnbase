context("OnDiskMSnExp class")

library(msdata)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))
inMem <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE)
onDisk <- readMSData2(files = mzf, msLevel. = 1, centroided. = TRUE)

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
multiMsInMem1 <- readMSData(files = f, msLevel. = 1, centroided. = TRUE)
multiMsInMem2 <- readMSData(files = f, msLevel. = 2, centroided. = TRUE)
multiMsOnDisk <- readMSData2(files = f, centroided. = TRUE)

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
    f <- msdata:::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    x <- readMSData2(f, verbose = FALSE)
    y <- readMSData(f, msLevel. = 2, centroided. = NA, verbose = FALSE)
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
    res2 <- MSnbase:::spectrapply(onDisk, ionCount)
    expect_identical(res1, res2)

    ## Plain function
    res1 <- lapply(spl, function(z) return(mean(mz(z))))
    res2 <- MSnbase:::spectrapply(onDisk, function(z) {
        return(mean(mz(z)))
    })
    expect_identical(res1, res2)

    ## Additional arguments.
    res1 <- lapply(spl, function(z, int) {
        return(mean(mz(z)[intensity(z) > int]))
    }, int = 30)
    res2 <- MSnbase:::spectrapply(onDisk, function(z, int) {
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
    sp1 <- multiMsInMem1[[13]]
    sp2 <- onDisk1[[13]]
    expect_identical(sp1, sp2)
    sp1 <- multiMsInMem1[[22]]
    sp2 <- onDisk1[[22]]
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
