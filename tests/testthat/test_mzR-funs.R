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