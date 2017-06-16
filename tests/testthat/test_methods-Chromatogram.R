test_that("Chromatogram accessors", {
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- sort(rnorm(100, mean = 300, sd = 3))
    chr <- Chromatogram(intensity = int, rt = rt)
    expect_equal(intensity(chr), int)
    expect_equal(rtime(chr), rt)
    chr@aggregationFun <- "sum"
    expect_true(validObject(chr))
    expect_equal(aggregationFun(chr), "sum")
    chr@aggregationFun <- "max"
    expect_true(validObject(chr))
    expect_equal(aggregationFun(chr), "max")
    chr@aggregationFun <- "mean"
    expect_true(validObject(chr))
    expect_equal(aggregationFun(chr), "mean")
    chr@aggregationFun <- "other"
    expect_error(validObject(chr))
    ## fromFile
    chr@aggregationFun <- "max"
    chr@fromFile <- 3L
    expect_true(validObject(chr))
    chr <- Chromatogram(intensity = int, rt = rt, fromFile = 2L)
    expect_true(validObject(chr))
    expect_equal(fromFile(chr), 2L)
    ## length
    expect_equal(length(chr), length(rt))
    ## as.data.frame
    df <- as.data.frame(chr)
    expect_equal(df, data.frame(rtime = rt, intensity = int))
    ## mz
    chr@mz <- c(1, 3)
    expect_equal(mz(chr), c(1, 3))
    chr <- Chromatogram(intensity = int, rt = rt, mz = 5)
    expect_equal(mz(chr), c(5, 5))
    expect_equal(mz(chr, filter = TRUE), c(NA_real_, NA_real_))
    chr@filterMz <- 4
    expect_equal(mz(chr, filter = TRUE), 4)
    chr <- Chromatogram(intensity = int, rt = rt, filterMz = 5)
    expect_equal(mz(chr, filter = TRUE), c(5, 5))
    ## precursorMz
    chr@precursorMz <- 123
    expect_equal(precursorMz(chr), 123)
    chr <- Chromatogram(intensity = int, rt = rt, precursorMz = 123)
    expect_equal(precursorMz(chr), c(123, 123))
    ## productMz
    chr@productMz <- 432
    expect_equal(productMz(chr), 432)
    chr <- Chromatogram(intensity = int, rt = rt, productMz = 432)
    expect_equal(productMz(chr), c(432, 432))
})

test_that("clean,Chromatogram works", {
    ints <- c(0, 3, 5, 0, 0, 0, 8, 9, 0, 9)
    rtime <- 1:length(ints)
    chr <- Chromatogram(rtime = rtime, intensity = ints)
    chr_2 <- clean(chr, all = TRUE)
    expect_equal(intensity(chr_2), c(3, 5, 8, 9, 9))
    expect_equal(rtime(chr_2), c(2, 3, 7, 8, 10))
    chr_2 <- clean(chr, all = FALSE)
    expect_equal(intensity(chr_2), c(0, 3, 5, 0, 0, 8, 9, 0, 9))
})
