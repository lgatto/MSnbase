context("MSmap class")

skip_on_travis()

library("rpx")
px1 <- PXDataset("PXD000001")
f <- pxget(px1, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")
ms <- mzR::openMSfile(f)
hd <- header(ms)

## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
    hd$retentionTime[ms1] / 60 < 35

test_that("MSmap accessors", {
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005)
    expect_null(show(M))

    expect_identical(mzRes(M), 0.005)
    expect_identical(dim(msMap(M)), c(75L, 401L))
    expect_identical(dim(M), c(75L, 401L))
    expect_identical(nrow(M), 75L)
    expect_identical(ncol(M), 401L)
    expect_identical(fileNames(M), fileName(M))

    sel <- hd$retentionTime/60 > 30 &
        hd$retentionTime/60 < 35 &
        hd$msLevel == 1
    expect_identical(hd$retentionTime[sel], rtime(M))
    expect_identical(msLevel(M), rep(1L, 75))
    expect_identical(mz(M), seq(521, 523, 0.005))

    Mt <- t(M)
    expect_identical(msMap(M), t(msMap(Mt)))
    expect_null(show(Mt))
})


test_that("MSmap from mzRpwiz and OnDiskMSnExp", {
    ## Use one of the sciex files.
    msn <- filterFile(sciex, 1)
    ## select slice
    rtsel <- which(rtime(msn) > 30 & rtime(msn) < 35)
    ## Map with mzRpwiz
    fh <- openMSfile(fileNames(msn))
    hd <- header(fh)
    M0 <- MSmap(fh, rtsel, 110, 112, .005, hd)
    mzR::close(fh)
    ## Map with OnDiskMSnExp
    M1 <- MSmap(msn, rtsel, 110, 112, .005)
    ## compare
    M0@call <- M1@call
    expect_equal(M0, M1, check.attributes = FALSE)
})

test_that("map data.frame", {
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005)
    mdf <- as(M, "data.frame")
    expect_equal(nrow(mdf), 401 * 75)
    expect_equal(colnames(mdf), c("intensity", "rt", "mz", "ms"))
    k <- sample(nrow(mdf), 1)
    i <- which(rtime(M)/60 == mdf[k, "rt"])
    j <- which(mz(M) == mdf[k, "mz"])
    expect_identical(msMap(M)[i, j], mdf[k, "intensity"])
})


## test_that("MSmap plotting", {
##     M <- MSmap(ms, ms1[rtsel], 521, 523, .005)
##     x <- plot3D(M)
##     expect_is(x, "trellis")
##     x <- plot(M, aspect = 1, allTicks = FALSE)
##     expect_is(x, "trellis")
## })
