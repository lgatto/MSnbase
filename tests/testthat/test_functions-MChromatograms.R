test_that("MChromatograms works", {
    chs <- new("MChromatograms")
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    expect_equal(unname(nrow(chs@phenoData)), 0)
    expect_equal(unname(nrow(chs@featureData)), 0)
    chs <- MChromatograms()
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    expect_equal(unname(nrow(chs@phenoData)), 0)
    expect_equal(unname(nrow(chs@featureData)), 0)
    ## Errors:
    chs@.Data <- matrix(1:4)
    expect_error(validObject(chs))
    expect_error(MChromatograms(1:4))
    ## Real data
    ch <- new("Chromatogram")
    ch_list <- list(ch, ch, ch, ch, ch, ch, ch, ch)
    chs <- MChromatograms(ch_list, nrow = 2)
    expect_equal(length(chs[1, 1]), 0)
    chs[2, 1] <- list(Chromatogram(1:10, 1:10))
    expect_equal(length(chs[2, 1]), 10)
    expect_equal(chs[, 1, drop = TRUE], c(`1` = ch, `2` = Chromatogram(1:10, 1:10)))
    expect_equal(unname(nrow(chs@phenoData)), ncol(chs))
    expect_equal(colnames(chs), rownames(pData(chs)))
    expect_equal(colnames(chs), as.character(1:ncol(chs)))
    expect_equal(rownames(chs), as.character(1:nrow(chs)))
    ## MChromatograms with a phenoData.
    pheno <- AnnotatedDataFrame(data.frame(idx = 1:4, name = letters[1:4]))
    chs <- MChromatograms(ch_list, nrow = 2, phenoData = pheno)
    expect_equal(chs@phenoData, as(pheno, "AnnotatedDataFrame"))
    expect_equal(colnames(chs), as.character(1:ncol(chs)))
    rownames(pheno) <- letters[1:4]
    chs <- MChromatograms(ch_list, nrow = 2, phenoData = pheno)
    expect_equal(chs@phenoData, as(pheno, "AnnotatedDataFrame"))
    expect_equal(colnames(chs), letters[1:4])
    ## MChromatograms with a featureData.
    fd <- data.frame(name = c("a", "b"), mzmin = 1:2)
    chs <- MChromatograms(ch_list, nrow = 2, featureData = fd)
    expect_equal(chs@featureData, AnnotatedDataFrame(fd))
    expect_equal(chs, MChromatograms(ch_list, nrow = 2,
                                    featureData = AnnotatedDataFrame(fd)))
    expect_equal(rownames(chs), as.character(1:nrow(chs)))
    rownames(fd) <- letters[1:2]
    chs <- MChromatograms(ch_list, nrow = 2, featureData = fd)
    expect_equal(rownames(chs), letters[1:2])
    expect_equal(featureNames(chs), letters[1:2])
    ## Error
    pheno <- AnnotatedDataFrame(data.frame(idx = 1:3, name = letters[1:3]))
    expect_error(MChromatograms(ch_list, nrow = 2, phenoData = pheno))
    expect_error(MChromatograms(ch_list, nrow = 2, featureData = 1:3))
    expect_error(MChromatograms(ch_list, nrow = 2,
                               featureData = data.frame(a = 1:5)))
})

test_that(".validMChromatograms works", {
    ch <- new("Chromatogram")
    ch_list <- list(ch, ch, ch, ch, ch, ch, ch, ch)
    chs <- MChromatograms(ch_list, nrow = 2)
    expect_true(MSnbase:::.validMChromatograms(chs))
    chs@phenoData <- as(AnnotatedDataFrame(), "AnnotatedDataFrame")
    expect_error(validObject(chs))
    chs <- MChromatograms(ch_list, nrow = 2)
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
    chrs <- MChromatograms(list(ch1, ch2, ch3, ch4, ch1, ch2), ncol = 2,
                          byrow = TRUE)
    MSnbase:::.plotChromatogramList(chrs[1, ])
    MSnbase:::.plotChromatogramList(chrs[1, ], main = "some title",
                                    col = c("red", "blue"))
    expect_error(MSnbase:::.plotChromatogramList(1:10))
})

test_that("c,MChromatograms, .bind_rows_chromatograms works", {
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))

    res <- MSnbase:::.bind_rows_chromatograms(chrs)
    expect_equal(res, chrs)
    res <- c(chrs)
    expect_equal(res, chrs)
    res <- MSnbase:::.bind_rows_chromatograms(chrs, chrs)
    expect_true(nrow(res) == 6)
    expect_equal(res[3, ], chrs[3, ])
    expect_equal(res[6, ], chrs[3, ])

    res_2 <- c(chrs, chrs)
    expect_equal(res, res_2)
    expect_true(nrow(res) == 6)
    expect_equal(res[3, ], chrs[3, ])
    expect_equal(res[6, ], chrs[3, ])

    res <- MSnbase:::.bind_rows_chromatograms(list(chrs, chrs))
    expect_equal(res, res_2)

    chrs <- MChromatograms(list(chr1, chr2, chr3), ncol = 3)

    res <- c(chrs, chrs)
    expect_equal(chrs[1, 2], res[1, 2])
    expect_equal(chrs[1, 3], res[2, 3])

    chrs_2 <- MChromatograms(list(chr1, chr2), ncol = 2)
    expect_error(c(chrs, chrs_2), "must match")
})
