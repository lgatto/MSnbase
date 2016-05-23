context("MSmap class")

library("AnnotationHub")
ah <- AnnotationHub()
ms <- ah[["AH49008"]]
hd <- header(ms)

test_that("", {
    ## a set of spectra of interest: MS1 spectra eluted
    ## between 30 and 35 minutes retention time
    ms1 <- which(hd$msLevel == 1)
    rtsel <- hd$retentionTime[ms1] / 60 > 30 &
        hd$retentionTime[ms1] / 60 < 35
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005)
    
    expect_identical(dim(msMap(M)), c(75L, 401L))
    expect_identical(dim(M), c(75L, 401L))
    expect_identical(nrow(M), 75L)
    expect_identical(ncol(M), 401L)
    expect_identical(fileNames(M), fileName(M))


})
