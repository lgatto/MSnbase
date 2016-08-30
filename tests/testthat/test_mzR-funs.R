test_that("list2Spectrum2 function", {
    data(itraqdata)
    sp <- itraqdata[[1]]
    l <- list(int = intensity(sp),
              mz = mz(sp))
    sp2 <- MSnbase:::list2Spectrum2(l)
    ## comparing intensities and mz only
    expect_identical(as(sp, "data.frame"),
                     as(sp2, "data.frame"))
    expect_identical(peaksCount(sp), peaksCount(sp2))
    expect_equal(tic(sp), tic(sp2))
    expect_identical(msLevel(sp), msLevel(sp2))
})

test_that("peaksAsList function", {
    f <- dir(system.file("extdata", package = "MSnbase"),
             full.names = TRUE, pattern = "mzXML")
    ms <- openMSfile(f)
    pl <- MSnbase:::peaksAsLists(ms)
    pli <- MSnbase:::peaksAsLists(ms, i = 1)
    plii <- MSnbase:::peaksAsLists(ms, i = 1:2)
    expect_identical(length(pl), length(ms))
    expect_length(pli, 1)
    expect_length(plii, 2)
    expect_identical(pli, pl[1])
    expect_identical(plii, pl[1:2])
})

test_that(".chomatogram function", {
    f <- dir(system.file("extdata", package = "MSnbase"),
             full.names = TRUE, pattern = "mzXML")
    ms <- openMSfile(f)
    hd <- header(ms)
    x <- MSnbase:::.chromatogram(hd, plot = FALSE)
    expect_is(x, "data.frame")
    expect_identical(nrow(x), length(ms))
    expect_identical(colnames(x), c("rt", "tic"))
    expect_identical(x$rt, hd$retentionTime)
    expect_equal(x$tic, 100 * hd$totIonCurrent/ max(hd$totIonCurrent))    
})


test_that("chomatogram methods", {
    f <- dir(system.file("microtofq", package = "msdata"),
             full.names = TRUE, pattern = "MM8.mzML")
    library("mzR")
    ms <- openMSfile(f)
    ch1 <- chromatogram(f, plot = FALSE)
    ch2 <- chromatogram(ms, plot = FALSE)
    expect_identical(ch1, ch2)
    hd <- header(ms)
    ch3 <- MSnbase:::.chromatogram(hd, plot = FALSE)
    expect_identical(ch1, ch3)
})

test_that("xic", {
    f <- dir(system.file("microtofq", package = "msdata"),
             full.names = TRUE, pattern = "MM8.mzML")
    ms <- openMSfile(f)

    expect_error(xic(ms))
    expect_warning(xicres <- xic(ms, mz = 636.925, width = 0.01, plot = FALSE))
    xicres1 <- xic(ms, mz = 636.925)
    xicres2 <- xic(f, mz = 636.925)
    expect_identical(xicres1, xicres2)
})
