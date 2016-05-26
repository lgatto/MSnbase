context("MSmap class")

library("AnnotationHub")
ah <- AnnotationHub()
ms <- ah[["AH49008"]]
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


test_that("MSmap plotting", {
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005)
    x <- plot3D(M)
    expect_is(x, "trellis")
    x <- plot(M, aspect = 1, allTicks = FALSE)
    expect_is(x, "trellis")
})
