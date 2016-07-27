############################################################
## Compare the "conventional" constructor with the C implementation
ints <- 1:1000
mzs <- 1:1000

test_that("Spectrum1 constructor", {
    system.time(
        Res1 <- new("Spectrum1", intensity = ints, mz = mzs, polarity = 1L,
                    rt = 12.4, fromFile = 3L, tic = 1234.3, centroided = TRUE)
    ) ## 0.009 sec
    ## The C constructor.
    system.time(
        Res2 <- MSnbase:::Spectrum1(intensity = ints, mz = mzs, polarity = 1L,
                                    rt = 12.4, fromFile = 3L, tic = 1234.3,
                                    centroided = TRUE)
    ) ## 0?
    expect_identical(Res1, Res2)

    ## Test exception, i.e. mz specified but not intensity or vice versa.
    expect_error(Test <- MSnbase:::Spectrum1(intensity = ints, polarity = 1L,
                                             rt = 12.4, fromFile = 3L,
                                             tic = 1234.3, centroided = TRUE))
    expect_error(Test <- MSnbase:::Spectrum1(mz = mzs, polarity = 1L, rt = 12.4,
                                             fromFile = 3L, tic = 1234.3,
                                             centroided = TRUE))
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
    system.time(
        sp1 <- new("Spectrum1", intensity = intVals, mz = mzVals,
                   polarity = 1L, fromFile = 1L, rt = 13.3, tic = 1234.3)
    ) ## 0.004
    ## unsorted C-constructor:
    system.time(
        sp2 <- MSnbase:::Spectrum1(intensity = intVals, mz = mzVals,
                                   polarity = 1L, fromFile = 1L, rt = 13.3,
                                   tic = 1234.3)
    )
    expect_identical(mz(sp1), sort(mz(sp2)))
    ## C-constructor with sorting:
    system.time(
        sp3 <- MSnbase:::Spectrum1_mz_sorted(intensity = intVals, mz = mzVals,
                                             polarity = 1L, fromFile = 1L,
                                             rt = 13.3, tic = 1234.3)
    ) ## 0.000
    expect_identical(mz(sp3), sort(mzVals))
    expect_identical(intensity(sp3), intVals[idx])
    expect_identical(sp3, sp1)
})

## Test the c-level multi-Spectrum1 constructor with M/Z ordering.
## o Ensure that ordering is as expected.
## o Check that the Spectrum values are as expected.
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
        spectL <- MSnbase:::Spectra1_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals)
    ) ## 0.003
    ## Check the TIC: should be the sum of intensities
    ticL <- lapply(intValsList, sum)
    expect_equal(unname(unlist(ticL)),
                 unname(unlist(lapply(spectL, tic))))
    ## The same with tic specified
    system.time(
        spectL <- MSnbase:::Spectra1_mz_sorted(rt = rts, acquisitionNum = acqN,
                                               scanIndex = acqN, mz = mzVals,
                                               intensity = intVals,
                                               fromFile = rep(1, 4),
                                               polarity = rep(1L, 4),
                                               nvalues = nvals,
                                               tic = rep(12, 4))
    ) ## 0.005
    expect_equal(rep(12, 4),
                 unname(unlist(lapply(spectL, tic))))

    ## Check if we've got the M/Z and intensity values correctly sorted.
    for (i in 1:length(idxList)) {
        expect_identical(mz(spectL[[i]]), mzValsList[[i]][idxList[[i]]])
        expect_identical(intensity(spectL[[i]]), intValsList[[i]][idxList[[i]]])
    }
})
