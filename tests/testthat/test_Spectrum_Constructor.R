############################################################
## Compare the "conventional" constructor with the C implementation
ints <- 1:1000
mzs <- 1:1000

test_that("Spectrum1 constructor", {
    system.time(
        Res1 <- new("Spectrum1", intensity=ints, mz=mzs, polarity=1L, rt=12.4, fromFile=3L,
                    tic=1234.3, centroided=TRUE)
    ) ## 0.009 sec
    ## The C constructor.
    system.time(
        Res2 <- MSnbase:::Spectrum1(intensity=ints, mz=mzs, polarity=1L, rt=12.4, fromFile=3L,
                                    tic=1234.3, centroided=TRUE)
    ) ## 0?
    expect_identical(Res1, Res2)
    ## Test exception, i.e. mz specified but not intensity or vice versa.
    expect_error(Test <- MSnbase:::Spectrum1(intensity=ints, polarity=1L, rt=12.4, fromFile=3L,
                                     tic=1234.3, centroided=TRUE))
    expect_error(Test <- MSnbase:::Spectrum1(mz=mzs, polarity=1L, rt=12.4, fromFile=3L,
                                     tic=1234.3, centroided=TRUE))
})

