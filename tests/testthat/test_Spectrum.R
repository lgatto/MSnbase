context("Spectrum class")

test_that("Spectrum validity", {
  expect_true(validObject(new("Spectrum1")))
  expect_true(validObject(new("Spectrum2")))
  expect_true(validObject(sp3 <- new("Spectrum2", msLevel = 3L)))
  expect_equal(msLevel(sp3), 3L)
})

test_that("Empty spectrum after trimming/filterrMz", {
  int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
           0, 1, 5, 10, 5, 1)
  sp <- new("Spectrum2",
            intensity = int,
            mz = 1:length(int))
  expect_true(validObject(sp))
  expect_warning(emptysp <- filterMz(sp, c(100, 110)))
  expect_true(validObject(emptysp))
  expect_true(isEmpty(emptysp))
})

test_that("Spectrum processing", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    sp <- new("Spectrum2",
              intensity = int,
              mz = 1:length(int))
    centroided(sp) <- FALSE
    ## removePeaks
    defaultT <- min(intensity(sp)[intensity(sp)>0])
    sp2a <- removePeaks(sp)
    sp2b <- removePeaks(sp, defaultT)
    sp2c <- removePeaks(sp, 3)
    expect_true(identical(sp2a, sp2b))
    expect_true(identical(sp, sp2b))
    expect_false(identical(sp, sp2c))
    expect_equal(ionCount(sp), sum(int))
    expect_equal(peaksCount(sp2c), peaksCount(sp))
    expect_equal(ionCount(sp), 55)
    expect_equal(ionCount(sp2c), 45)
    ## clean
    sp3 <- clean(sp)
    expect_equal(ionCount(sp), ionCount(sp3))
    expect_equal(peaksCount(sp), length(int))
    expect_equal(peaksCount(sp3), 23)
    ##trimMz
    sp4 <- trimMz(sp, c(10, 20))
    expect_equal(intensity(sp4), int[10:20])
    expect_equal(mz(sp4), 10:20)
    expect_equal(peaksCount(sp4), length(10:20))
    expect_equal(ionCount(sp4), sum(int[10:20]))
})

test_that("Spectrum normalisation", {
    s1 <- new("Spectrum1", mz = 1:5, intensity = 1:5)
    s2 <- new("Spectrum2", mz = 1:5, intensity = 1:5,
              precursorIntensity = 10)
    ## Spectrum1
    ## max is default
    expect_equal(intensity(normalize(s1)), (1:5) / 5)
    expect_equal(intensity(normalise(s1)), (1:5) / 5)
    expect_equal(intensity(normalize(s1, method = "max")), (1:5) / 5)
    expect_equal(intensity(normalize(s1, method = "sum")), (1:5) / 15)
    expect_error(normalize(s1, method = "precursor"), "'arg' should be one of")
    ## Spectrum2
    ## max is default
    expect_equal(intensity(normalize(s2)), (1:5) / 5)
    expect_equal(intensity(normalise(s2)), (1:5) / 5)
    expect_equal(intensity(normalize(s2, method = "max")), (1:5) / 5)
    expect_equal(intensity(normalize(s2, method = "sum")), (1:5) / 15)
    expect_equal(intensity(normalize(s2, method = "precursor")), (1:5) / 10)
    expect_equal(intensity(normalize(s2, method = "precursor",
                                     precursorIntensity = 20)), (1:5) / 20)
})

test_that("Noise estimation", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = c(1:3, 2:1))
    s2 <- new("Spectrum2", mz = 3, intensity = 3, centroided = TRUE)
    e <- matrix(NA, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(estimateNoise(new("Spectrum2")), "spectrum is empty")
    expect_warning(estimateNoise(s2), "only supported for profile spectra")
    expect_equal(suppressWarnings(estimateNoise(new("Spectrum2"))), e)
    expect_equal(suppressWarnings(estimateNoise(s2)), e)
    centroided(s1) <- FALSE
    expect_equal(estimateNoise(s1),
                 cbind(mz = 1:5, intensity = mad(intensity(s1))))
})

test_that("Peak picking", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = c(1:3, 2:1))
    s2 <- new("Spectrum2", mz = 3, intensity = 3, centroided = TRUE)
    expect_warning(pickPeaks(new("Spectrum2")), "spectrum is empty")
    expect_equal(suppressWarnings(pickPeaks(new("Spectrum2"))),
                 new("Spectrum2"))
    expect_equal(suppressWarnings(pickPeaks(s2)), s2)
    centroided(s1) <- FALSE
    expect_equal(pickPeaks(s1), s2)
})

test_that("Spectrum smoothing", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = c(1:3, 2:1))
    s2 <- new("Spectrum2", mz = 1:5, intensity = c(2, 2, 2+1/3, 2, 2))
    expect_warning(smooth(new("Spectrum2")), "spectrum is empty")
    expect_equal(smooth(s1, method = "MovingAverag", halfWindowSize = 1),
                 s2)
})

test_that("Spectrum quantification", {
    ## dummy Spectrum
    int <- c(0, 2, 3, 1, 0)
    mz <- c(114.11,
            114.12,
            114.13,
            114.14,
            114.145)
    sp <- new("Spectrum2",
              intensity = int,
              mz = mz,
              centroided = FALSE)
    expect_true(validObject(sp))
    expect_equal(MSnbase:::getCurveWidth(sp, iTRAQ4[1]),
                 list(lwr = 1, upr = 5))
    expect_equal(as.numeric(quantify(sp, "sum", iTRAQ4[1])$peakQuant), 6)
    expect_equal(as.numeric(quantify(sp, "max", iTRAQ4[1])$peakQuant), 3)
    expect_that(as.numeric(quantify(sp, "trap", iTRAQ4[1])$peakQuant),
                equals((0.01 * 2) / 2 +
                       (0.01 * 2)     +
                       (0.01 * 1) / 2 +
                       0.01 * 1       +
                       (0.01 * 2) / 2 +
                       (0.01 * 0.5) / 2))
    ## print("Warnings expected because there is not data for
    ## iTRAQ4[2].") -- not since v1.1.2
    expect_true(as.logical(is.na(quantify(sp, "sum", iTRAQ4[2])$peakQuant)))
    ## expect_warning(quantify(sp,"sum",iTRAQ4[2])$peakQuant)
})

test_that("Spectrum strict quantification", {
    ## dummy Spectrum
    int <- c(0, 1, 1, 3, 1, 1, 0)
    mz. <- c(113.9,
             114.0,
             114.05,
             114.1,
             114.15,
             114.2,
             114.25)
    sp <- new("Spectrum2",
              intensity = int,
              mz = mz.,
              centroided = FALSE)
    expect_true(validObject(sp))
    expect_equivalent(
        quantify(sp, "trap", iTRAQ4[1], strict = FALSE)$peakQuant,
        (mz.[2] - mz.[1]) * (int[2] - int[1]) / 2 +
        (mz.[3] - mz.[2]) * int[3] +
        (mz.[4] - mz.[3]) * int[3] +
        (mz.[4] - mz.[3]) * (int[4] - int[3]) / 2 +
        (mz.[5] - mz.[4]) * int[5] +
        (mz.[5] - mz.[4]) * (int[4] - int[5]) / 2 +
        (mz.[6] - mz.[5]) * int[6] +
        (mz.[7] - mz.[6]) * (int[6] - int[7]) / 2)
    ## changing width to keep calculation below correct, since
    ## reporter ions mz changed in commit
    ## c82c82bd20af5840375abca0f7f41f7f36e8e4ef
    iTRAQ4@width <- 0.065
    expect_equivalent(
        quantify(sp, "trap", iTRAQ4[1], strict = TRUE)$peakQuant ,
        (mz.[4] - mz.[3]) * int[3] +
        (mz.[4] - mz.[3]) * (int[4] - int[3]) / 2 +
        (mz.[5] - mz.[4]) * int[5] +
        (mz.[5] - mz.[4]) * (int[4] - int[5]) / 2)
})

## test_that("breaks_Spectrum", {
##     s1 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
##     ## issue 191
##     expect_equal(MSnbase:::breaks_Spectrum(s1), 1:5)
##     expect_equal(MSnbase:::breaks_Spectrum(s1, breaks = 1:2), c(1, 2, 5))
##     expect_equal(MSnbase:::breaks_Spectrum(s1, binSize = 2), c(1, 3, 6))
## })

test_that(".fix_breaks works as breaks_Spectra", {
    s1 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    s2 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    brks <- seq(floor(min(c(mz(s1), mz(s1)))),
                ceiling(max(c(mz(s1), mz(s1)))), by = 1)
    expect_equal(brks, 1:4)
    expect_equal(MSnbase:::.fix_breaks(brks, c(1, 4)), 1:5)
    brks <- seq(floor(min(c(mz(s1), mz(s2)))),
                ceiling(max(c(mz(s1), mz(s2)))), by = 1)
    expect_equal(brks, 1:5)
    ## issue 190
    expect_equal(MSnbase:::.fix_breaks(brks, c(1, 5)), 1:6)
    brks <- seq(floor(min(c(mz(s1), mz(s2)))),
                ceiling(max(c(mz(s1), mz(s2)))), by = 2)
    expect_equal(brks, c(1, 3, 5))
    expect_equal(MSnbase:::.fix_breaks(brks, c(1, 6)), c(1, 3, 5, 7))

    s3 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    s4 <- new("Spectrum2", mz = 11:15, intensity = 1:5)
    brks <- seq(floor(min(c(mz(s3), mz(s4)))),
                ceiling(max(c(mz(s3), mz(s4)))), by = 1)
    expect_equal(brks, 1:15)
    expect_equal(MSnbase:::.fix_breaks(brks, c(1, 15)), 1:16)
    brks <- seq(floor(min(c(mz(s3), mz(s4)))),
                ceiling(max(c(mz(s3), mz(s4)))), by = 2)
    expect_equal(brks, seq(1, 15, 2))
    expect_equal(MSnbase:::.fix_breaks(brks, c(1, 15)), seq(1, 17, by=2))
})

test_that("bin_Spectrum", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    s2 <- new("Spectrum2", mz = 1:5 + 0.1, intensity = 1:5)
    r1 <- new("Spectrum2", mz = 1:5 + 0.5, intensity = 1:5, tic = 15)
    r2 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(3, 7, 5), tic = 15)
    r3 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5), tic = 10)
    r31 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5), tic = 10)
    r4 <- new("Spectrum2", mz = c(1, 3, 5), intensity = c(1, 5, 9), tic = 15)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 1), r1)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 2), r2)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 2, fun = mean), r3)
    expect_equal(MSnbase:::bin_Spectrum(s1, breaks = seq(0, 7, by = 2)), r4)
    expect_equal(MSnbase:::bin_Spectrum(s2, binSize = 1), r1)
    expect_equal(MSnbase:::bin_Spectrum(s2, binSize = 2, fun = mean), r31)
    expect_equal(MSnbase:::bin_Spectrum(s2, breaks = seq(0, 7, by = 2)), r4)
})

test_that("bin_Spectrum - bug fix #ecaaa324505b17ee8c4855806f7e37f14f1b27b8", {
    s <- new("Spectrum2", mz = c(1:7, 55, 78, 100), intensity = 1:10)
    s2 <- bin(s)
    expect_equal(mz(s2), c(seq(1.5, 100.5, 1)))
    ires <- rep(0, peaksCount(s2))
    ires[peaksCount(s2)] <- intensity(s)[peaksCount(s)]
    ires[1:7] <- 1:7
    ires[55] <- 8
    ires[78] <- 9
    expect_equal(intensity(s2), ires)
})

test_that("bin_Spectra", {
    # issue 190
    s1 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    s2 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    r1 <- new("Spectrum2", mz = 1:5 + 0.5, intensity = c(1:4, 0))
    r2 <- new("Spectrum2", mz = 1:5 + 0.5, intensity = 1:5)
    r3 <- new("Spectrum2", mz = 1:4 + 0.5, intensity = 1:4)
    expect_equal(MSnbase:::bin_Spectra(s1, s2), list(r1, r2))
    expect_equal(MSnbase:::bin_Spectra(s1, s1), list(r3, r3))
})

test_that("removePeaks profile vs centroided", {
    int <- c(2, 0, 0, 0, 1, 5, 1, 0, 0, 1, 3, 1, 0, 0, 1, 4, 2, 1)
    sp1 <- new("Spectrum2",
               intensity = int,
               centroided = FALSE,
               mz = 1:length(int))
    res1 <- c(0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    expect_identical(intensity(removePeaks(sp1, 4)), res1)

    res2 <- int <- c(104, 57, 32, 33, 118, 76, 38, 39, 52, 140, 52,
                     88, 394, 71, 408, 94, 2032)
    sp2 <- new("Spectrum2",
               intensity = int,
               centroided = FALSE,
               mz = seq_len(length(int)))
    expect_identical(intensity(removePeaks(sp2, 500)),
                     intensity(sp2))
    res2[res2 < 500] <- 0

    expect_identical(intensity(removePeaks(sp2, 500)),
                     intensity(sp2))
    centroided(sp2) <- TRUE
    expect_identical(intensity(removePeaks(sp2, 500)),
                     res2)
})

test_that("empty spectrum", {
    s <- new("Spectrum2")
    expect_true(isEmpty(s))
    t <- removePeaks(s, 10)
    expect_true(all.equal(s, t))
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6,
             2, 1, 0, 1, 2, 0, 0, 1, 5, 10, 5, 1)
    sp <- new("Spectrum2",
              intensity = int,
              mz = 1:length(int))
    expect_false(isEmpty(sp))
})

test_that("show MS1 spectrum", {
    f <- dir(system.file("threonine", package = "msdata"),
             full.names = TRUE)
    x <- readMSData(f, msLevel = 1)
    expect_null(show(x[[1]]))
})

test_that(".spectrum_header works", {
    mzf <- mzR::openMSfile(fileNames(tmt_erwinia_on_disk))
    hdr <- header(mzf)
    mzR::close(mzf)
    sp_1 <- tmt_erwinia_on_disk[[1]]
    sp_2 <- tmt_erwinia_on_disk[[2]]

    hdr_1 <- .spectrum_header(sp_1)
    expect_equal(unname(hdr_1["acquisitionNum"]), hdr$acquisitionNum[1])
    expect_equal(hdr$collisionEnergy[1], unname(hdr_1["collisionEnergy"]))
    expect_equal(hdr$highMZ[1], unname(hdr_1["highMZ"]))
    expect_equal(hdr$ionisationEnergy[1], unname(hdr_1["ionisationEnergy"]))
    expect_equal(hdr$lowMZ[1], unname(hdr_1["lowMZ"]))
    expect_equal(hdr$mergedResultEndScanNum[1],
                 unname(hdr_1["mergedResultEndScanNum"]))
    expect_equal(hdr$mergedResultScanNum[1],
                 unname(hdr_1["mergedResultScanNum"]))
    expect_equal(hdr$mergedResultStartScanNum[1],
                 unname(hdr_1["mergedResultStartScanNum"]))
    expect_equal(hdr$mergedScan[1], unname(hdr_1["mergedScan"]))
    expect_equal(hdr$msLevel[1], unname(hdr_1["msLevel"]))
    expect_equal(hdr$peaksCount[1], unname(hdr_1["peaksCount"]))
    expect_equal(hdr$polarity[1], unname(hdr_1["polarity"]))
    expect_equal(hdr$precursorCharge[1], unname(hdr_1["precursorCharge"]))
    expect_equal(hdr$precursorIntensity[1], unname(hdr_1["precursorIntensity"]))
    expect_equal(hdr$precursorMZ[1], unname(hdr_1["precursorMZ"]))
    expect_equal(hdr$precursorScanNum[1], unname(hdr_1["precursorScanNum"]))
    expect_equal(hdr$retentionTime[1], unname(hdr_1["retentionTime"]))
    ## Failing: base peak most likely because the data was filtered,
    ## injectionTime because it's not stored in a Spectrum object
    ## expect_equal(hdr$basePeakIntensity[1], unname(hdr_1["basePeakIntensity"]))
    ## expect_equal(hdr$basePeakMZ[1], unname(hdr_1["basePeakMZ"]))
    ## expect_equal(hdr$injectionTime[1], unname(hdr_1["injectionTime"]))

    hdr_2 <- .spectrum_header(sp_2)
    expect_equal(unname(hdr_2["acquisitionNum"]), hdr$acquisitionNum[2])
    expect_equal(hdr$collisionEnergy[2], unname(hdr_2["collisionEnergy"]))
    expect_equal(hdr$highMZ[2], unname(hdr_2["highMZ"]))
    expect_equal(hdr$ionisationEnergy[2], unname(hdr_2["ionisationEnergy"]))
    expect_equal(hdr$lowMZ[2], unname(hdr_2["lowMZ"]))
    expect_equal(hdr$mergedResultEndScanNum[2],
                 unname(hdr_2["mergedResultEndScanNum"]))
    expect_equal(hdr$mergedResultScanNum[2],
                 unname(hdr_2["mergedResultScanNum"]))
    expect_equal(hdr$mergedResultStartScanNum[2],
                 unname(hdr_2["mergedResultStartScanNum"]))
    expect_equal(hdr$mergedScan[2], unname(hdr_2["mergedScan"]))
    expect_equal(hdr$msLevel[2], unname(hdr_2["msLevel"]))
    expect_equal(hdr$peaksCount[2], unname(hdr_2["peaksCount"]))
    expect_equal(hdr$polarity[2], unname(hdr_2["polarity"]))
    expect_equal(hdr$precursorCharge[2], unname(hdr_2["precursorCharge"]))
    expect_equal(hdr$precursorIntensity[2], unname(hdr_2["precursorIntensity"]))
    expect_equal(hdr$precursorMZ[2], unname(hdr_2["precursorMZ"]))
    expect_equal(hdr$precursorScanNum[2], unname(hdr_2["precursorScanNum"]))
    expect_equal(hdr$retentionTime[2], unname(hdr_2["retentionTime"]))
    ## Failing: base peak most likely because the data was filtered,
    ## injectionTime because it's not stored in a Spectrum object
    ## expect_equal(hdr$basePeakIntensity[1], unname(hdr_1["basePeakIntensity"]))
    ## expect_equal(hdr$basePeakMZ[1], unname(hdr_1["basePeakMZ"]))
    ## expect_equal(hdr$injectionTime[1], unname(hdr_1["injectionTime"]))
})

test_that("kNeighbors works", {
    ## Test the m/z refining method for peak picking/centroiding.
    ints <- c(3, 4, 5, 7, 3, 4, 2, 8, 5, 6, 8, 8.1, 4, 5, 6, 3)
    mzs <- 1:length(ints) + rnorm(length(ints), mean = 0, sd = 0.1)
    plot(mzs, ints, type = "h")
    pk_pos <- c(4, 8, 12)

    res <- kNeighbors(mzs, ints, peakIdx = pk_pos, k = 1)
    points(res[, 1], res[, 2], type = "h", col = "blue")
    expect_equal(unname(res[1, 1]), weighted.mean(mzs[3:5], ints[3:5]))
    expect_equal(unname(res[2, 1]), weighted.mean(mzs[7:9], ints[7:9]))
    expect_equal(unname(res[3, 1]), weighted.mean(mzs[11:13], ints[11:13]))

    res <- kNeighbors(mzs, ints, peakIdx = pk_pos, k = 2)
    points(res[, 1], res[, 2], type = "h", col = "green")
    expect_equal(unname(res[1, 1]), weighted.mean(mzs[2:6], ints[2:6]))
    expect_equal(unname(res[2, 1]), weighted.mean(mzs[6:10], ints[6:10]))
    expect_equal(unname(res[3, 1]), weighted.mean(mzs[10:14], ints[10:14]))
    
    expect_error(kNeighbors(mz = 3, ints))
})

test_that("descendPeak works", {
    ints <- c(2, 3, 1, 2, 1, 0, 1, 2, 0, 1, 0, 2, 3, 2, 1, 2, 5, 8, 7, 6,
              5, 4, 3, 2, 1, 0, 1, 1, 4)
    mzs <- 1:length(ints) + rnorm(length(ints), mean = 0, sd = 0.1)
    plot(mzs, ints, type = "h")
    pk_pos <- c(13, 18)

    res <- descendPeak(mzs, ints, pk_pos, signalPercentage = 0)
    points(res[, 1], res[, 2], type = "h", col = "blue")
    expect_equal(unname(res[1, 1]), weighted.mean(mzs[11:15], ints[11:15]))
    expect_equal(unname(res[2, 1]), weighted.mean(mzs[15:26], ints[15:26]))
    
    res <- descendPeak(mzs, ints, pk_pos, signalPercentage = 0,
                       stopAtTwo = TRUE)
    points(res[, 1], res[, 2], type = "h", col = "green")
    expect_equal(unname(res[1, 1]), weighted.mean(mzs[6:15], ints[6:15]))
    expect_equal(unname(res[2, 1]), weighted.mean(mzs[15:26], ints[15:26]))

    ## With signalPercentage
    res <- descendPeak(mzs, ints, pk_pos, signalPercentage = 50,
                       stopAtTwo = TRUE)
    points(res[, 1], res[, 2], type = "h", col = "orange")
    idx <- 6:15
    idx <- idx[ints[idx] > ints[13]/2]
    expect_equal(unname(res[1, 1]), weighted.mean(mzs[idx], ints[idx]))
    idx <- 15:26
    idx <- idx[ints[idx] > ints[18]/2]
    expect_equal(unname(res[2, 1]), weighted.mean(mzs[idx], ints[idx]))
})

test_that("pickPeaks,Spectrum works with refineMz", {
    ## Get one spectrum from the tmt
    spctr <- tmt_erwinia_in_mem_ms1[[1]]
    centroided(spctr) <- FALSE

    mzr <- c(530.9, 531.2)
    plot(mz(filterMz(spctr, mz = mzr)), intensity(filterMz(spctr, mz = mzr)),
         type = "h")
    ## plain pickPeaks
    spctr_pks <- pickPeaks(spctr)
    points(mz(filterMz(spctr_pks, mz = mzr)),
           intensity(filterMz(spctr_pks, mz = mzr)),
           type = "p", col = "blue")
    ## Now the same but using a refineMz method.
    spctr_kn <- pickPeaks(spctr, refineMz = "kNeighbors", k = 1)
    points(mz(filterMz(spctr_kn, mz = mzr)),
           intensity(filterMz(spctr_kn, mz = mzr)),
           type = "p", col = "red")
    ## Now the same but using a refineMz method.
    spctr_kn <- pickPeaks(spctr, refineMz = "kNeighbors", k = 2)
    points(mz(filterMz(spctr_kn, mz = mzr)),
           intensity(filterMz(spctr_kn, mz = mzr)),
           type = "p", col = "green")
    spctr_kn <- pickPeaks(spctr, refineMz = "descendPeak",
                          signalPercentage = 45)
    points(mz(filterMz(spctr_kn, mz = mzr)),
           intensity(filterMz(spctr_kn, mz = mzr)),
           type = "p", col = "red")
    spctr_kn <- pickPeaks(spctr, refineMz = "descendPeak",
                          signalPercentage = 10, stopAtTwo = TRUE)
    points(mz(filterMz(spctr_kn, mz = mzr)),
           intensity(filterMz(spctr_kn, mz = mzr)),
           type = "p", col = "orange")
        
    ## Check if we can call method and refineMz and pass arguments to both
    spctr_kn <- pickPeaks(spctr, refineMz = "kNeighbors", k = 1,
                          method = "SuperSmoother", span = 0.9)
    points(mz(filterMz(spctr_kn, mz = mzr)),
           intensity(filterMz(spctr_kn, mz = mzr)),
           type = "p", col = "red")
    
    ## Check errors
    expect_error(pickPeaks(spctr, refineMz = "some_method"))
    expect_error(pickPeaks(spctr, not_sup = TRUE, method = "SuperSmoother"))
})

test_that(".combineMovingWindow works for Spectrum", {
    ## on a list of spectra.
    spcts <- spectra(tmt_erwinia_in_mem_ms1)
    s_comb <- .combineMovingWindow(spcts)
    expect_equal(length(spcts), length(s_comb))
    expect_equal(unname(lapply(spcts, rtime)), lapply(s_comb, rtime))
    expect_equal(unname(lapply(spcts, msLevel)), lapply(s_comb, msLevel))

    ## Check the first.
    vals_exp <- do.call(rbind, lapply(spcts[1:2], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[1]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[1]]), vals_exp$i)
    
    ## Check the second.
    vals_exp <- do.call(rbind, lapply(spcts[1:3], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[2]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[2]]), vals_exp$i)

    ## With halfWindowSize 4L
    s_comb <- .combineMovingWindow(spcts, halfWindowSize = 4L)
    expect_equal(length(spcts), length(s_comb))
    expect_equal(unname(lapply(spcts, rtime)), lapply(s_comb, rtime))
    expect_equal(unname(lapply(spcts, msLevel)), lapply(s_comb, msLevel))

    ## Check the first.
    vals_exp <- do.call(rbind, lapply(spcts[1:5], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[1]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[1]]), vals_exp$i)
    
    ## Check the fifth
    vals_exp <- do.call(rbind, lapply(spcts[1:9], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[5]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[5]]), vals_exp$i)
})

test_that(".estimate_mz_scattering works", {
    set.seed(123)
    mzs <- seq(1, 20, 0.1)
    all_mz <- c(mzs + rnorm(length(mzs), sd = 0.01),
                mzs + rnorm(length(mzs), sd = 0.005),
                mzs + rnorm(length(mzs), sd = 0.02))
    res <- .estimate_mz_scattering(sort(all_mz))
    expect_true(length(res) == 1)
    expect_true(res < 0.051)
    expect_error(.estimate_mz_scattering(mzs))

    all_mz <- c(mzs + rnorm(length(mzs), sd = 0.01),
                mzs + rnorm(length(mzs), sd = 0.005),
                mzs + rnorm(length(mzs), sd = 0.06))
    res <- .estimate_mz_scattering(sort(all_mz))
    expect_true(res < 0.08)
    expect_true(length(res) == 1)
})

test_that(".group_mz_values works", {
    set.seed(123)
    mzs <- seq(1, 20, 0.1)
    all_mz <- sort(c(mzs + rnorm(length(mzs), sd = 0.001),
                     mzs + rnorm(length(mzs), sd = 0.005),
                     mzs + rnorm(length(mzs), sd = 0.002)))
    res <- MSnbase:::.group_mz_values(all_mz)
    expect_true(length(res) == length(all_mz))
    ## Expect groups of 3 each.
    expect_true(all(table(res) == 3))
    
    ## Remove one from the 2nd group.
    res <- MSnbase:::.group_mz_values(all_mz[-5])
    expect_true(sum(res == 2) == 2)
})

test_that("combineSpectra works", {
    set.seed(123)
    mzs <- seq(1, 20, 0.1)
    mzs_2 <- c(mzs, 20.1)
    ints1 <- abs(rnorm(length(mzs), 10))
    ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
    ints2 <- c(abs(rnorm(length(mzs), 10)), 4)
    ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
    ints3 <- abs(rnorm(length(mzs), 10))
    ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)

    ## Create the spectra
    sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
               intensity = ints1, rt = 1)
    sp2 <- new("Spectrum1", mz = mzs_2 + rnorm(length(mzs_2), sd = 0.01),
               intensity = ints2, rt = 2)
    sp3 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.008),
               intensity = ints3, rt = 3)
    sp4 <- new("Spectrum2", mz = mzs + rnorm(length(mzs), sd = 0.3),
               intensity = ints2, rt = 4)
    expect_error(combineSpectra(list(sp1, sp2, sp3, sp4)))

    res <- combineSpectra(list(sp1, sp2, sp3), timeDomain = TRUE)
    expect_equal(length(mz(res)), length(mz(sp2)))
    expect_equal(rtime(res), rtime(sp2))

    res <- combineSpectra(list(sp2, sp1), timeDomain = FALSE)
    expect_equal(length(mz(res)), length(mz(sp1)))
    expect_equal(rtime(res), rtime(sp1))

    sp4 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.3),
               intensity = ints2, rt = 4)
    ## randon noise larger than resolution.
    expect_error(res <- combineSpectra(list(sp1, sp3, sp4)))

    res <- combineSpectra(list(sp1, sp2, sp3), main = 1, timeDomain = TRUE)
    expect_equal(rtime(res), rtime(sp1))
    expect_equal(length(mz(res)), length(mz(sp1)))

    res <- combineSpectra(list(sp1, sp2, sp3), main = 3, timeDomain = TRUE)
    expect_equal(rtime(res), rtime(sp3))
    expect_equal(length(mz(res)), length(mz(sp3)))
    
    res <- combineSpectra(list(sp1, sp1), intensityFun = sum, timeDomain = TRUE)
    expect_equal(mz(res), mz(sp1))
    expect_equal(intensity(res), intensity(sp1) * 2)

    ## Use character mzFun:
    expect_error(combineSpectra(list(sp1, sp2, sp3), mzFun = "meani"))
    res <- combineSpectra(list(sp1, sp2, sp3), mzFun = base::mean)
    res2 <- combineSpectra(list(sp1, sp2, sp3), mzFun = "weighted.mean")
    expect_equal(intensity(res), intensity(res2))
    expect_false(all(mz(res) == mz(res2)))

    ## Use real data.
    od1 <- filterFile(sciex, 1)
    lst <- spectra(od1[3:5])

    res <- combineSpectra(lst, timeDomain = TRUE)
    res_2 <- combineSpectra(lst, timeDomain = FALSE)

    expect_equal(mz(res), mz(res_2))
    expect_equal(intensity(res), intensity(res_2))
    ## with (wrongly) pre-calculated mzd
    mzd <- MSnbase:::.estimate_mz_scattering(sort(unlist(lapply(lst, mz))))
    expect_error(combineSpectra(lst, timeDomain = TRUE, mzd = mzd))
    res_3 <- combineSpectra(lst, timeDomain = FALSE, mzd = mzd)
    
    expect_equal(mz(res), mz(res_3))
    expect_equal(intensity(res), intensity(res_3))
})

test_that(".estimate_mz_resolution, estimateMzResolution,Spectrum works", {
    set.seed(123)
    mzs <- seq(1, 2000, 0.1)
    mzs <- mzs + rnorm(length(mzs), sd = 0.005)
    res <- .estimate_mz_resolution(mzs)
    ## expect_true(res - 0.1 < 0.005)
    expect_true(res - 0.1 < 0.008)

    res1 <- estimateMzResolution(tmt_erwinia_in_mem_ms1[[1]])
    res2 <- estimateMzResolution(tmt_erwinia_in_mem_ms1[[2]])
    expect_true(res1 != res2)
})

test_that(".findPeakValley works", {
    vals <- c(3, 5, 6, 7, 8, 9, 5, 4, 2, 1, 5, 7, 4, 1)
    expect_equal(.findPeakValley(6:20, vals), 10)
    expect_equal(.findPeakValley(6:1, vals), NA)
    expect_equal(.findPeakValley(12:14, vals), NA)
    expect_equal(.findPeakValley(12:1, vals), 10)
})

test_that(".density works", {
    set.seed(123)
    xs <- rnorm(300, 2, 45)
    res <- .density(xs)
    res_2 <- density(xs, n = 512L)
    expect_equal(res$x, res_2$x)
    expect_equal(res$y, res_2$y)
})

