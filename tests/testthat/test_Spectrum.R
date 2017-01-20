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


test_that("bin_Spectrum", {
    s1 <- new("Spectrum2", mz = 1:5, intensity = 1:5)
    s2 <- new("Spectrum2", mz = (1:5)+0.1, intensity = 1:5)
    r1 <- new("Spectrum2", mz = c(1.5, 2.5, 3.5, 4.5, 5.5), intensity = 1:5, tic = 15)
    r2 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(3, 7, 5), tic = 15)
    r3 <- new("Spectrum2", mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5), tic = 10)
    r4 <- new("Spectrum2", mz = c(1, 3, 5), intensity = c(1, 5, 9), tic = 15)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 1), r1)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 2), r2)
    expect_equal(MSnbase:::bin_Spectrum(s1, binSize = 2, fun = mean), r3)
    expect_equal(MSnbase:::bin_Spectrum(s1, breaks = seq(0, 7, by = 2)), r4)
    expect_equal(MSnbase:::bin_Spectrum(s2, binSize = 1), r1)
    expect_equal(MSnbase:::bin_Spectrum(s2, binSize = 2), r2)
    expect_equal(MSnbase:::bin_Spectrum(s2, binSize = 2, fun = mean), r3)
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
