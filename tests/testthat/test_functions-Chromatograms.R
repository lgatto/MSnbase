test_that("Chromatograms works", {
    chs <- new("Chromatograms")
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    expect_equal(unname(nrow(chs@phenoData)), 0)
    expect_equal(unname(nrow(chs@featureData)), 0)
    chs <- Chromatograms()
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    expect_equal(unname(nrow(chs@phenoData)), 0)
    expect_equal(unname(nrow(chs@featureData)), 0)
    ## Errors:
    chs@.Data <- matrix(1:4)
    expect_error(validObject(chs))
    expect_error(Chromatograms(1:4))
    ## Real data
    ch <- new("Chromatogram")
    ch_list <- list(ch, ch, ch, ch, ch, ch, ch, ch)
    chs <- Chromatograms(ch_list, nrow = 2)
    expect_equal(length(chs[1, 1]), 0)
    chs[2, 1] <- list(Chromatogram(1:10, 1:10))
    expect_equal(length(chs[2, 1]), 10)
    expect_equal(chs[, 1, drop = TRUE], c(`1` = ch, `2` = Chromatogram(1:10, 1:10)))
    expect_equal(unname(nrow(chs@phenoData)), ncol(chs))
    expect_equal(colnames(chs), rownames(pData(chs)))
    expect_equal(colnames(chs), as.character(1:ncol(chs)))
    expect_equal(rownames(chs), as.character(1:nrow(chs)))
    ## Chromatograms with a phenoData.
    pheno <- AnnotatedDataFrame(data.frame(idx = 1:4, name = letters[1:4]))
    chs <- Chromatograms(ch_list, nrow = 2, phenoData = pheno)
    expect_equal(chs@phenoData, as(pheno, "NAnnotatedDataFrame"))
    expect_equal(colnames(chs), as.character(1:ncol(chs)))
    rownames(pheno) <- letters[1:4]
    chs <- Chromatograms(ch_list, nrow = 2, phenoData = pheno)
    expect_equal(chs@phenoData, as(pheno, "NAnnotatedDataFrame"))
    expect_equal(colnames(chs), letters[1:4])
    ## Chromatograms with a featureData.
    fd <- data.frame(name = c("a", "b"), mzmin = 1:2)
    chs <- Chromatograms(ch_list, nrow = 2, featureData = fd)
    expect_equal(chs@featureData, AnnotatedDataFrame(fd))
    expect_equal(chs, Chromatograms(ch_list, nrow = 2,
                                    featureData = AnnotatedDataFrame(fd)))
    expect_equal(rownames(chs), as.character(1:nrow(chs)))
    rownames(fd) <- letters[1:2]
    chs <- Chromatograms(ch_list, nrow = 2, featureData = fd)
    expect_equal(rownames(chs), letters[1:2])
    expect_equal(featureNames(chs), letters[1:2])
    ## Error
    pheno <- AnnotatedDataFrame(data.frame(idx = 1:3, name = letters[1:3]))
    expect_error(Chromatograms(ch_list, nrow = 2, phenoData = pheno))
    expect_error(Chromatograms(ch_list, nrow = 2, featureData = 1:3))
    expect_error(Chromatograms(ch_list, nrow = 2,
                               featureData = data.frame(a = 1:5)))
})

test_that(".validChromatograms works", {
    ch <- new("Chromatogram")
    ch_list <- list(ch, ch, ch, ch, ch, ch, ch, ch)
    chs <- Chromatograms(ch_list, nrow = 2)
    expect_true(MSnbase:::.validChromatograms(chs))
    chs@phenoData <- as(AnnotatedDataFrame(), "NAnnotatedDataFrame")
    expect_error(validObject(chs))
    chs <- Chromatograms(ch_list, nrow = 2)
    chs@featureData <- AnnotatedDataFrame()
    expect_error(validObject(chs))
})

test_that(".plotChromatogramList works", {
    ints <- abs(rnorm(123, mean = 200, sd = 19))
    ch1 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
    ints <- abs(rnorm(122, mean = 300, sd = 35))
    ch2 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
    ints <- abs(rnorm(124, mean = 214, sd = 49))
    ch3 <- Chromatogram(rtime = seq_along(ints) + 300, intensity = ints,
                        mz = 403)
    ints <- abs(rnorm(123, mean = 530, sd = 89))
    ch4 <- Chromatogram(rtime = seq_along(ints) + 300, intensity = ints,
                        mz = 403)
    chrs <- Chromatograms(list(ch1, ch2, ch3, ch4, ch1, ch2), ncol = 2,
                          byrow = TRUE)
    MSnbase:::.plotChromatogramList(chrs[1, ])
    MSnbase:::.plotChromatogramList(chrs[1, ], main = "some title",
                                    col = c("red", "blue"))
    expect_error(MSnbase:::.plotChromatogramList(1:10))
})
