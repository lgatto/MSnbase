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

test_that("breaks_Spectrum", {
    s1 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    ## issue 191
    expect_equal(MSnbase:::breaks_Spectrum(s1), 1:5)
    expect_equal(MSnbase:::breaks_Spectrum(s1, breaks = 1:2), c(1, 2, 5))
    expect_equal(MSnbase:::breaks_Spectrum(s1, binSize = 2), c(1, 3, 6))
})

test_that("breaks_Spectra", {
    s1 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    s2 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    expect_equal(MSnbase:::breaks_Spectra(s1, s1), 1:5)
    expect_equal(MSnbase:::breaks_Spectra(s2, s2), 1:6)
    ## issue 190
    expect_equal(MSnbase:::breaks_Spectra(s1, s2), 1:6)
    expect_equal(MSnbase:::breaks_Spectra(s1, s2, binSize = 2), c(1, 3, 5, 7))

    s3 <- new("Spectrum2", mz = 1:4, intensity = 1:4)
    s4 <- new("Spectrum2", mz = 11:15, intensity = 1:5)
    expect_equal(MSnbase:::breaks_Spectra(s3, s4), 1:16)
    expect_equal(MSnbase:::breaks_Spectra(s3, s4, binSize = 2), seq(1, 17, by=2))
})

test_that("bin_Spectrum", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    s2 <- new("Spectrum2", mz = 1:5 + 0.1, intensity = 1:5)
    r1 <- new("Spectrum2", mz = 1:5 + 0.5, intensity = 1:5, tic = 15)
    r2 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(3, 7, 5), tic = 15)
    r3 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5), tic = 10)
    r31 <- new("Spectrum2", mz = c(2, 4, 6.5), intensity = c(1.5, 3.5, 5), tic = 10)
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
