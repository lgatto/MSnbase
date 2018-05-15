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

    ## msLevel
    chr@msLevel <- 2L
    expect_equal(msLevel(chr), 2L)
    expect_true(validObject(chr))
    chr@msLevel <- 1:4
    expect_equal(msLevel(chr), 1:4)
    expect_true(validObject(chr))

    chr <- Chromatogram(intensity = int, rt = rt, msLevel = 4)
    expect_equal(msLevel(chr), 4L)
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
    
    chr <- Chromatogram(
        rtime = 1:12,
        intensity = c(0, 0, 20, 0, 0, 0, 123, 124343, 3432, 0, 0, 0))
    chr_clnd <- clean(chr)
    expect_equal(rtime(chr_clnd), c(2, 3, 4, 6, 7, 8, 9,10))
    expect_equal(intensity(chr_clnd), c(0, 20, 0, 0, 123, 124343, 3432, 0))
    chr_clnd <- clean(chr, all = TRUE)
    expect_true(length(chr_clnd) == 4)
    expect_equal(rtime(chr_clnd), c(3, 7, 8, 9))

    ## With NA
    chr <- Chromatogram(
        rtime = 1:12,
        intensity = c(0, NA, 20, 0, 0, 0, 123, 124343, 3432, 0, 0, 0))
    chr_clnd <- clean(chr)
    expect_equal(intensity(chr_clnd), c(NA, 20, 0, 0, 123, 124343, 3432, 0))
    expect_equal(rtime(chr_clnd), c(2, 3, 4, 6, 7, 8, 9, 10))
    chr_clnd <- clean(chr, na.rm = TRUE)
    expect_equal(intensity(chr_clnd), c(0, 20, 0, 0, 123, 124343, 3432, 0))
    expect_equal(rtime(chr_clnd), c(1, 3, 4, 6, 7, 8, 9, 10))
    chr <- Chromatogram(
        rtime = 1:12,
        intensity = c(NA, NA, 20, NA, NA, NA, 123, 124343, 3432, NA, NA, NA))
    chr_clnd <- clean(chr)
    expect_equal(intensity(chr_clnd), c(NA, 20, NA, NA, 123, 124343, 3432, NA))
    expect_equal(rtime(chr_clnd), c(2, 3, 4, 6, 7, 8, 9, 10))
    chr_clnd <- clean(chr, na.rm = TRUE)
    expect_equal(intensity(chr_clnd), c(20, 123, 124343, 3432))
    expect_equal(rtime(chr_clnd), c(3, 7, 8, 9))
})

test_that("filterRt,Chromatogram works", {
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    chr <- Chromatogram(intensity = int, rtime = sort(rt))

    chr_2 <- filterRt(chr, rt = c(200, 300))
    expect_true(all(rtime(chr_2) >= 200))
    expect_true(all(rtime(chr_2) <= 300))
    ints <- intensity(chr_2)
    expect_equal(ints, intensity(chr)[rtime(chr) >= 200 & rtime(chr) <= 300])

    ## No rt
    expect_equal(chr, filterRt(chr))

    ## Outside range
    chr_2 <- filterRt(chr, rt = c(400, 500))
    expect_true(length(chr_2) == 0)
    expect_equal(intensity(chr_2), numeric())
    expect_equal(rtime(chr_2), numeric())
})

test_that("isEmpty,Chromatogram and plot,Chromatogram work", {
    chr <- Chromatogram()
    expect_true(isEmpty(chr))
    expect_warning(plot(chr))
    
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    chr <- Chromatogram(intensity = int, rtime = sort(rt))
    expect_true(!isEmpty(chr))
    plot(chr)
    
    chr <- Chromatogram(intensity = rep_len(NA_real_, length(rt)),
                        rtime = sort(rt))
    expect_true(isEmpty(chr))
    expect_warning(plot(chr))
})

test_that(".bin_Chromatogram and bin,Chromatogram work", {
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- seq(2, 200, length.out = 100)
    chr <- Chromatogram(intensity = int, rtime = rt)

    chrb <- MSnbase:::.bin_Chromatogram(chr, binSize = 4)
    expect_equal(rtime(chrb), seq(4, 200, by = 4))
    expect_equal(bin(chr, binSize = 4), chrb)
    expect_error(.bin_Chromatogram(chr, fun = "bls"))

    ## Provide breaks that are larger than rtime.
    brks <- seq(2, 300, by = 4)
    res <- .bin_Chromatogram(chr, breaks = brks)
    expect_equal(length(rtime(res)), (length(brks) - 1))
    expect_equal(intensity(res)[1:length(chrb)], intensity(chrb))
})
