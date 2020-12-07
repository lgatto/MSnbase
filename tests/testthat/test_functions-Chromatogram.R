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

test_that(".align_chromatogram_approx works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8),
                         intensity = c(5, 9, 3, 1, 4, 3, 6, 9))
    chr2 <- Chromatogram(rtime = c(3, 4, 6), intensity = c(3, 1, 3))
    res <- .align_chromatogram_approx(chr2, chr1)
    expect_equal(length(chr1), length(res))
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, NA, 3, 1, 2, 3, NA, NA))

    res <- .align_chromatogram_approx(chr1, chr2)
    expect_equal(length(chr2), length(res))
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(3, 1, 3))

    ## Not perfectly matching rtimes:
    chr1 <- Chromatogram(rtime = c(1.1, 2.1, 3.1, 4.1, 5.1),
                         intensity = c(1, 2, 3, 2, 1))
    chr2 <- Chromatogram(rtime = c(2, 3), intensity = c(3, 5))
    res <- .align_chromatogram_approx(chr2, chr1)
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, 3.2, NA, NA, NA))

    res <- .align_chromatogram_approx(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(1.9, 2.9))
})

test_that(".align_chromatogram_closest works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8),
                         intensity = c(5, 9, 3, 1, 4, 3, 6, 9))
    chr2 <- Chromatogram(rtime = c(3, 4, 6), intensity = c(3, 1, 3))
    res <- .align_chromatogram_closest(chr2, chr1)
    expect_equal(length(chr1), length(res))
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, NA, 3, 1, NA, 3, NA, NA))

    res <- .align_chromatogram_closest(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(3, 1, 3))
    res_none <- .align_chromatogram_closest(chr1, chr2, tolerance = 0)
    expect_equal(res, res_none)

    ## Not perfectly matching rtimes:
    chr1 <- Chromatogram(rtime = c(1.1, 2.1, 3.1, 4.1, 5.1),
                         intensity = c(1, 2, 3, 2, 1))
    chr2 <- Chromatogram(rtime = c(2, 3), intensity = c(3, 5))
    res <- .align_chromatogram_closest(chr2, chr1)
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, 3, 5, NA, NA))
    res <- .align_chromatogram_closest(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(2, 3))

    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8),
                         intensity = c(5, 9, 3, 1, 4, 3, 6, 9))
    chr2 <- Chromatogram(rtime = c(3, 4, 6), intensity = c(3, 1, 3))
    res <- .align_chromatogram_closest(chr2, chr1, tolerance = 0)
    expect_equal(length(chr1), length(res))
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, NA, 3, 1, NA, 3, NA, NA))

    ## Not perfectly matching rtimes:
    chr1 <- Chromatogram(rtime = c(1.1, 2.1, 3.1, 4.1, 5.1),
                         intensity = c(1, 2, 3, 2, 1))
    chr2 <- Chromatogram(rtime = c(2, 3), intensity = c(3, 5))
    res <- .align_chromatogram_closest(chr2, chr1, tolerance = 0)
    expect_equal(rtime(res), rtime(chr1))
    expect_true(all(is.na(intensity(res))))
})

test_that(".align_chromatogram works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
    chr2 <- Chromatogram(rtime = c(2.5, 3.42, 4.5, 5.43, 6.5),
                         intensity = c(5, 12, 15, 11, 5))
    res <- .align_chromatogram(chr1, chr2)
    expect_equal(res, .align_chromatogram_closest(chr1, chr2))
    res <- .align_chromatogram(chr1, chr2, method = "approx")
    expect_equal(res, .align_chromatogram_approx(chr1, chr2))

    expect_error(.align_chromatogram(chr1, chr2, method = "other"),
                 "should be one of")
})
