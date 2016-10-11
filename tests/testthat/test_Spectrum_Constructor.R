############################################################
## Compare the "conventional" constructor with the C implementation
ints <- as.numeric(1:1000)
mzs <- as.numeric(1:1000)

############################################################
## Spectrum1.
test_that("Spectrum1 constructor", {
    Res1 <- new("Spectrum1", intensity = ints, mz = mzs, polarity = 1L,
                rt = 12.4, fromFile = 3L, tic = 1234.3, centroided = TRUE)
    Res2 <- MSnbase:::Spectrum1(intensity = ints, mz = mzs, polarity = 1L,
                                rt = 12.4, fromFile = 3L, tic = 1234.3,
                                centroided = TRUE)
    expect_identical(Res1, Res2)
    ## Test exception, i.e. mz specified but not intensity or vice versa.
    expect_error(Test <- MSnbase:::Spectrum1(intensity = ints, polarity = 1L,
                                             rt = 12.4, fromFile = 3L,
                                             tic = 1234.3, centroided = TRUE))
    expect_error(Test <- MSnbase:::Spectrum1(mz = mzs, polarity = 1L, rt = 12.4,
                                             fromFile = 3L, tic = 1234.3,
                                             centroided = TRUE))
    expect_identical(classVersion(Res1)["Spectrum1"],
                     new("Versions",
                         Spectrum1 = MSnbase:::getClassVersionString("Spectrum1")))
    expect_identical(classVersion(Res1)["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))

})

test_that("M/Z sorted Spectrum1 constructor", {
    mzVals <- abs(rnorm(3500, mean = 100, sd = 10))
    intVals <- abs(rnorm(3500, mean = 10, sd = 5))

    ##sorted <- MSnbase:::sortNumeric(mzVals)
    ##idx <- MSnbase:::orderNumeric(mzVals)
    ## expect_identical(idx, order(mzVals))
    sorted <- sort(mzVals)
    idx <- order(mzVals)
    ## R constructor:
    sp1 <- new("Spectrum1", intensity = intVals, mz = mzVals,
               polarity = 1L, fromFile = 1L, rt = 13.3, tic = 1234.3)
    ## unsorted C-constructor:
    sp2 <- MSnbase:::Spectrum1(intensity = intVals, mz = mzVals,
                                   polarity = 1L, fromFile = 1L, rt = 13.3,
                                   tic = 1234.3)
    expect_identical(mz(sp1), sort(mz(sp2)))
    ## C-constructor with sorting:
    sp3 <- MSnbase:::Spectrum1_mz_sorted(intensity = intVals, mz = mzVals,
                                         polarity = 1L, fromFile = 1L,
                                         rt = 13.3, tic = 1234.3)
    expect_identical(mz(sp3), sort(mzVals))
    expect_identical(intensity(sp3), intVals[idx])
    expect_identical(sp3, sp1)
    expect_identical(classVersion(sp3)["Spectrum1"],
                     new("Versions",
                         Spectrum1 = MSnbase:::getClassVersionString("Spectrum1")))
    expect_identical(classVersion(sp3)["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))

})

## Test the c-level multi-Spectrum1 constructor with M/Z ordering.
## o Ensure that ordering is as expected.
## o Check that the Spectrum values are as expected.
## o Check that classVersions are properly set.
test_that("C-level multi-Spectrum1 constructor with M/Z ordering", {
    ## Use Spectra1_mz_sorted constructor.
    mzVals <- abs(rnorm(20000, mean = 10, sd = 10))
    intVals <- abs(rnorm(20000, mean = 1000, sd = 1000))
    rts <- c(1.3, 1.4, 1.5, 1.6)
    acqN <- 1:4
    nvals <- rep(5000, 4)

    mzValsList <- split(mzVals, f = rep(1:4, each = 5000))
    intValsList <- split(intVals, f = rep(1:4, each = 5000))
    idxList <- lapply(mzValsList, order)

    system.time(
        ## Switch on gctorture to force potential memory mapping problems.
        ## gctorture(on = TRUE)
        spectL <- MSnbase:::Spectra1_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals)
        ## gctorture(on = FALSE)
    ) ## 0.003
    expect_true(all(unlist(lapply(spectL, validObject))))
    ## Check the TIC: should be the sum of intensities
    ticL <- lapply(intValsList, sum)
    expect_equal(unname(unlist(ticL)),
                 unname(unlist(lapply(spectL, tic))))
    ## Check the class version for one of the spectra:
    expect_identical(classVersion(spectL[[3]])["Spectrum1"],
                     new("Versions",
                         Spectrum1 = MSnbase:::getClassVersionString("Spectrum1")))
    expect_identical(classVersion(spectL[[1]])["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))
    ## The same with tic specified
    system.time(
        ##gctorture(on = TRUE)
        spectL <- MSnbase:::Spectra1_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals,
                                               tic = rep(12, 4))
        ## gctorture(on = FALSE)
    ) ## 0.005
    expect_true(all(unlist(lapply(spectL, validObject))))
    expect_equal(rep(12, 4),
                 unname(unlist(lapply(spectL, tic))))
    expect_identical(classVersion(spectL[[3]])["Spectrum1"],
                     new("Versions",
                         Spectrum1 = MSnbase:::getClassVersionString("Spectrum1")))
    expect_identical(classVersion(spectL[[1]])["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))
    ## Check if we've got the M/Z and intensity values correctly sorted.
    for (i in 1:length(idxList)) {
        expect_identical(mz(spectL[[i]]), mzValsList[[i]][idxList[[i]]])
        expect_identical(intensity(spectL[[i]]), intValsList[[i]][idxList[[i]]])
    }
})


############################################################
## Spectrum2.

## M/Z sorted Spectrum2 constructor.
test_that("M/Z sorted Spectrum2 constructor", {
    mzVals <- abs(rnorm(3500, mean = 100, sd = 10))
    intVals <- abs(rnorm(3500, mean = 10, sd = 5))

    sorted <- sort(mzVals)
    idx <- order(mzVals)
    ## R constructor:
    system.time(
        sp1 <- new("Spectrum2", intensity = intVals, mz = mzVals,
                   polarity = -1L, fromFile = 1L, rt = 13.3, tic = 1234.3)
    ) ## 0.004
    expect_identical(mz(sp1), sort(mzVals))
    expect_identical(intensity(sp1), intVals[idx])
    ## Check other slot values...
    expect_identical(sp1@polarity, -1L)
    expect_identical(sp1@fromFile, 1L)
    expect_identical(sp1@rt, 13.3)
    expect_identical(sp1@tic, 1234.3)
    ## Check class versions
    expect_identical(classVersion(sp1)["Spectrum2"],
                     new("Versions",
                         Spectrum2 = MSnbase:::getClassVersionString("Spectrum2")))
    expect_identical(classVersion(sp1)["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))
    ## C-constructor with sorting:
    sp2 <- MSnbase:::Spectrum2_mz_sorted(intensity = intVals, mz = mzVals,
                                         polarity = -1L, fromFile = 1L,
                                         rt = 13.3, tic = 1234.3)
    expect_identical(mz(sp2), sort(mzVals))
    expect_identical(intensity(sp2), intVals[idx])
    ## Check class versions
    expect_identical(classVersion(sp2)["Spectrum2"],
                     new("Versions",
                         Spectrum2 = MSnbase:::getClassVersionString("Spectrum2")))
    expect_identical(classVersion(sp2)["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))

    ## Calculate tic within:
    sp1 <- new("Spectrum2", intensity = intVals, mz = mzVals)
    expect_equal(sp1@tic, sum(intVals))
    ## Test some exceptions...
    ## o Pass only intensity or mz
    expect_error(new("Spectrum2", intensity = intVals, polarity = 1L))
    expect_error(new("Spectrum2", mz = muVals, polarity = 1L))
})

## Test the c-level multi-Spectrum2 constructor with M/Z ordering.
## o Ensure that ordering is as expected.
## o Check that the Spectrum values are as expected.
test_that("C-level multi-Spectrum2 constructor with M/Z ordering", {
    ## Use Spectra2_mz_sorted constructor.
    mzVals <- abs(rnorm(20000, mean = 10, sd = 10))
    intVals <- abs(rnorm(20000, mean = 1000, sd = 1000))
    rts <- c(1.3, 1.4, 1.5, 1.6)
    acqN <- 1:4
    nvals <- rep(5000, 4)

    mzValsList <- split(mzVals, f = rep(1:4, each = 5000))
    intValsList <- split(intVals, f = rep(1:4, each = 5000))
    idxList <- lapply(mzValsList, order)

    system.time(
        ## gctorture(on = TRUE)
        spectL <- MSnbase:::Spectra2_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals)
        ## gctorture(on = FALSE)
    ) ## 0.003
    expect_true(all(unlist(lapply(spectL, validObject))))
    ## Check the TIC: should be the sum of intensities
    ticL <- lapply(intValsList, sum)
    expect_equal(unname(unlist(ticL)),
                 unname(unlist(lapply(spectL, tic))))
    ## Check class versions
    expect_identical(classVersion(spectL[[1]])["Spectrum2"],
                     new("Versions",
                         Spectrum2 = MSnbase:::getClassVersionString("Spectrum2")))
    expect_identical(classVersion(spectL[[3]])["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))
    ## The same with tic specified
    system.time(
        ## gctorture(on = TRUE)
        spectL <- MSnbase:::Spectra2_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals,
                                               tic = rep(12, 4))
        ## gctorture(on = FALSE)
    ) ## 0.005
    expect_true(all(unlist(lapply(spectL, validObject))))
    expect_equal(rep(12, 4),
                 unname(unlist(lapply(spectL, tic))))
    ## Check class versions
    expect_identical(classVersion(spectL[[1]])["Spectrum2"],
                     new("Versions",
                         Spectrum2 = MSnbase:::getClassVersionString("Spectrum2")))
    expect_identical(classVersion(spectL[[3]])["Spectrum"],
                     new("Versions",
                         Spectrum = MSnbase:::getClassVersionString("Spectrum")))

    ## Check if we've got the M/Z and intensity values correctly sorted.
    for (i in 1:length(idxList)) {
        expect_identical(mz(spectL[[i]]), mzValsList[[i]][idxList[[i]]])
        expect_identical(intensity(spectL[[i]]), intValsList[[i]][idxList[[i]]])
    }
})
