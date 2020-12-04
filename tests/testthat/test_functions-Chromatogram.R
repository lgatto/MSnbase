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

test_that(".normalize_chromatogram works", {
    chr <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7),
                        intensity = c(NA_real_, 13, 16, 22, 34, 15, 6))
    res <- .normalize_chromatogram(chr)
    expect_true(max(intensity(res), na.rm = TRUE) == 1)
    expect_true(is.na(intensity(res)[1]))
})

test_that(".filter_intensity_chromatogram works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
    res <- .filter_intensity_chromatogram(chr1, intensity = 4)
    expect_true(all(intensity(res) > 4))

    res <- .filter_intensity_chromatogram(
        chr1, intensity = function(x) x@intensity > max(x@intensity) / 2)
    expect_true(all(intensity(res) > max(intensity(res) / 2)))

    expect_error(
        .filter_intensity_chromatogram(chr1, intensity = function(x) TRUE),
        "expected result")
})
