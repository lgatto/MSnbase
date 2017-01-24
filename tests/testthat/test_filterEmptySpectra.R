test_that("filterEmptySpectra", {
    require("msdata")
    mzf <- system.file("microtofq/MM14.mzML", package = "msdata")
    
    ## ondisk
    od <- readMSData2(files = mzf, msLevel. = 1, centroided. = TRUE)
    od2 <- filterMz(od, mz = c(1005, 1006))
    expect_identical(length(od), length(od2))
    expect_identical(length(filterEmptySpectra(od2)), 0L)
    
    ## inmem
    im <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE)
    im2 <- filterMz(im, mz = c(1005, 1006))
    expect_identical(length(im), length(im2))
    expect_identical(length(filterEmptySpectra(im2)), 0L)
})
