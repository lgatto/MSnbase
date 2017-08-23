test_that("Chromatogram construction works", {
    ## Constructing the object
    ch <- new("Chromatogram")
    ch@mz <- 3
    expect_error(validObject(ch))
    ch@mz <- c(1, 3)
    ch@precursorMz <- 4
    expect_error(validObject(ch))
    ch@precursorMz <- c(4, 4)
    ch@productMz <- 5
    expect_error(validObject(ch))
    ##
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    ## check exceptions:
    expect_error(Chromatogram(intensity = int))
    chr <- Chromatogram()
    chr@rtime <- rt
    expect_error(validObject(chr))
    chr@intensity <- int
    expect_error(validObject(chr))
    ## xcms issue #145: values are ordered based on rtime
    idx <- order(rt)
    chr@rtime <- rt[idx]
    chr@intensity <- int[idx]
    expect_true(validObject(chr))
    chr_2 <- Chromatogram(intensity = int, rtime = rt)
    expect_equal(chr_2, chr)
    expect_equal(rtime(chr), sort(rt))
    expect_equal(intensity(chr), int[order(rt)])
    expect_error(Chromatogram(aggregationFun = "other"))
    chr@aggregationFun <- "max"
    expect_true(validObject(chr))

    chr_2 <- Chromatogram(intensity = int, rtime = rt, msLevel = 2L)
    expect_equal(chr_2@msLevel, 2L)
})

test_that(".plotChromatogram works", {
    chr <- Chromatogram(rtime = 1:20, 1:20)
    MSnbase:::.plotChromatogram(chr, main = "dummy")
    plot(chr)
})
