test_that("[,MChromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    ## with a non-empty pData.
    chrs2 <- chrs
    pData(chrs2) <- data.frame(id = 1:2)
    fData(chrs2) <- data.frame(ion = c("a", "b"))
    chrs2[, 1]

    ## o Subset using indices
    expect_true(is(chrs[1, 1], "Chromatogram"))
    expect_equal(chrs[1, 2], ch2)
    ##   extract a row
    expect_equal(chrs[1, , drop = TRUE], list(`1` = ch, `2` = ch2))
    expect_equal(chrs[1, , drop = FALSE], MChromatograms(list(ch, ch2), nrow = 1))
    ##   Test the default
    expect_equal(chrs[1, ], MChromatograms(list(ch, ch2), nrow = 1))
    ##   extract a column
    expect_equal(chrs[, 2, drop = TRUE], list(`1` = ch2, `2` = ch3))
    res <- chrs[, 2, drop = FALSE]
    res_exp <- MChromatograms(list(ch2, ch3), ncol = 1,
                             dimnames = list(NULL, "2"))
    ## Have to re-place the rownames of pheno data othewise we compare numeric
    ## against character
    rownames(pData(res)) <- rownames(pData(res))
    expect_equal(res, res_exp)
    ##   Repeat with colnames:
    colnames(chrs) <- c("a", "b")
    expect_true(is(chrs[1, 1], "Chromatogram"))
    expect_equal(chrs[1, 2], ch2)
    ##   extract a row
    expect_equal(chrs[1, , drop = TRUE], list(a = ch, b = ch2))
    res_exp <- MChromatograms(list(ch, ch2), nrow = 1)
    colnames(res_exp) <- c("a", "b")
    expect_equal(chrs[1, , drop = FALSE], res_exp)
    ##   Test the default
    expect_equal(chrs[1, ], res_exp)
    ##   extract a column
    expect_equal(chrs[, 2, drop = TRUE], list(`1` = ch2, `2` = ch3))
    res_exp <- MChromatograms(list(ch2, ch3), ncol = 1)
    colnames(res_exp) <- "b"
    expect_equal(chrs[, 2, drop = FALSE], res_exp)
    ## Check also the featureData
    res <- chrs[2, ]
    expect_equal(rownames(res), "2")
    expect_equal(featureNames(res), "2")

    ## o Subset using logical
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    expect_true(is(chrs[c(TRUE, FALSE), c(TRUE, FALSE)], "Chromatogram"))
    expect_equal(chrs[c(TRUE, FALSE), c(FALSE, TRUE)], ch2)
    ##   extract a row
    expect_equal(chrs[c(TRUE, FALSE), , drop = TRUE], list(`1` = ch, `2` = ch2))
    expect_equal(chrs[c(TRUE, FALSE), , drop = FALSE],
                 MChromatograms(list(ch, ch2), nrow = 1))
    expect_equal(chrs[c(TRUE, FALSE), ],
                 MChromatograms(list(ch, ch2), nrow = 1))
    ##   extract a column
    expect_equal(chrs[, c(FALSE, TRUE), drop = TRUE], list(`1` = ch2, `2` = ch3))
    res <- chrs[, c(FALSE, TRUE), drop = FALSE]
    rownames(pData(res)) <- rownames(pData(res))
    expect_equal(res, MChromatograms(list(ch2, ch3), ncol = 1,
                                    dimnames = list(NULL, "2")))
    ##   Repeat with colnames
    colnames(chrs) <- c("a", "b")
    expect_equal(chrs[c(TRUE, FALSE), , drop = TRUE], list(a = ch, b = ch2))
    res_exp <- MChromatograms(list(ch, ch2), nrow = 1)
    colnames(res_exp) <- c("a", "b")
    expect_equal(chrs[c(TRUE, FALSE), , drop = FALSE], res_exp)
    expect_equal(chrs[c(TRUE, FALSE), ], res_exp)
    ## extract a column
    expect_equal(chrs[, c(FALSE, TRUE), drop = TRUE], list(`1` = ch2, `2` = ch3))
    res_exp <- MChromatograms(list(ch2, ch3), ncol = 1)
    colnames(res_exp) <- "b"
    expect_equal(chrs[, c(FALSE, TRUE)], res_exp)
    expect_equal(chrs[, c(FALSE, TRUE)], res_exp)

    ## Subset using names
    expect_equal(chrs[, "a", drop = TRUE], list(`1` = ch, `2` = ch1))
    res_exp <- MChromatograms(list(ch, ch1), ncol = 1)
    colnames(res_exp) <- "a"
    expect_equal(chrs[, "a", drop = FALSE], res_exp)
    expect_equal(chrs[, "a"], res_exp)

    ## Check phenoData while subsetting.
    pd <- data.frame(name = letters[1:2], idx = 1:2)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    res <- chrs[, 2]
    pd_exp <- droplevels(pd[2, ])
    expect_equal(pData(res), pd_exp)
    rownames(pd) <- c("g", "h")
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    res <- chrs[, 2]
    pd_exp <- droplevels(pd[2, ])
    expect_equal(pData(res), pd_exp)

    ## Check featureData while subsetting
    fd <- data.frame(a = c("first", "second"), mz = c(2, 4))
    rownames(fd) <- fd$a
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd),
                          featureData = fd)
    expect_equal(rownames(chrs), rownames(fd))
    expect_equal(fData(chrs[, 1]), fd)
    expect_equal(fData(chrs[2, ]), droplevels(fd[2, ]))
})

test_that("[<-,MChromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    colnames(chrs) <- c("a", "b")

    ints <- abs(rnorm(94, sd = 200))
    ch4 <- Chromatogram(rtime = 1:length(ints), ints)

    ## errors
    expect_error(chrs[1:2, 1:2] <- list(ch4, ch4, ch4, ch4))
    expect_error(chrs["z", ] <- list(ch4, ch4))

    ## Single element.
    chrs[1, 2] <- ch4
    expect_equal(chrs[1, 2], ch4)
    chrs[, 2] <- list(ch2, ch3)
    expect_equal(chrs[, 2, drop = TRUE], list(`1` = ch2, `2` = ch3))
    chrs[2, 1] <- list(ch4)
    expect_equal(chrs[2, 1], ch4)

    chrs[, "a"] <- list(ch2, ch1)
    expect_equal(chrs[, 1, drop = TRUE], list(`1` = ch2, `2` = ch1))
    expect_error(chrs[, 1] <- list(ch, ch2, ch3))

    chrs[, c(TRUE, FALSE)] <- list(ch4, ch4)
    expect_equal(chrs[, 1, drop = TRUE], list(`1` = ch4, `2` = ch4))
})

test_that("plot,MChromatograms works", {
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
    plot(chrs)
    plot(chrs[1, , drop = FALSE])
    plot(chrs[1, 1, drop = FALSE])
    plot(chrs[1, ])
    plot(chrs[1, 1])
    plot(chrs[, 2])
})

test_that("colnames<-, sampleNames, sampleNames<-,MChromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    expect_equal(colnames(chrs), as.character(1:ncol(chrs)))
    expect_equal(sampleNames(chrs), as.character(1:ncol(chrs)))
    colnames(chrs) <- letters[1:ncol(chrs)]
    expect_equal(colnames(chrs), letters[1:ncol(chrs)])
    expect_equal(rownames(pData(chrs)), letters[1:ncol(chrs)])
    expect_equal(sampleNames(chrs), letters[1:ncol(chrs)])

    sampleNames(chrs) <- c("b", "z")
    expect_equal(colnames(chrs), c("b", "z"))
    ## Error
    expect_error(colnames(chrs) <- 1:4)
})

test_that("phenoData,pData,pData<-,MChromatograms works", {
    ## Check if we can access the phenoData.
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    pd_exp <- annotatedDataFrameFrom(matrix(ncol = 2, nrow = 2), byrow = FALSE)
    rownames(pData(pd_exp)) <- NULL
    pd_exp <- as(pd_exp, "AnnotatedDataFrame")
    expect_equal(phenoData(chrs), pd_exp)

    ## Check error when assigning a phenoData with different names
    pd <- data.frame(name = letters[1:2], idx = 1:2)
    rownames(pd) <- letters[1:2]
    expect_error(pData(chrs) <- pd)

    pd <- data.frame(name = letters[1:2], idx = 1:2)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    chrs_2 <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                            phenoData = pd)
    expect_equal(chrs, chrs_2)
    expect_equal(phenoData(chrs), as(AnnotatedDataFrame(pd), "AnnotatedDataFrame"))

    ## pData.
    expect_equal(pData(chrs), pd)

    pd_2 <- cbind(pd, other = 1:2)
    pData(chrs) <- pd_2
    expect_equal(pData(chrs), pd_2)

    rownames(pd_2) <- c("g", "h")
    expect_error(pData(chrs) <- pd_2)
    colnames(chrs) <- c("g", "h")
    pData(chrs) <- pd_2
    expect_equal(pData(chrs), pd_2)
    expect_equal(colnames(chrs), rownames(pd_2))

    ## $
    expect_equal(chrs$name, letters[1:2])
    expect_equal(chrs$idx, 1:2)

    chrs$idx <- c(2, 1)
    expect_equal(chrs$idx, c(2, 1))

    chrs$new_variable <- c("it", "works")
    expect_equal(chrs$new_variable, c("it", "works"))

    expect_error(chrs$new_variable <- 1:4)
    chrs$new_variable <- 1
})

test_that("rownames<-, featureNames, featureNames<-,MChromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    expect_equal(rownames(chrs), as.character(1:nrow(chrs)))
    expect_equal(featureNames(chrs), as.character(1:nrow(chrs)))
    rownames(chrs) <- letters[1:nrow(chrs)]
    expect_true(validObject(chrs))
    expect_equal(rownames(chrs), letters[1:nrow(chrs)])
    expect_equal(featureNames(chrs), letters[1:nrow(chrs)])
    expect_error(rownames(chrs) <- letters[1:20])
    featureNames(chrs) <- c("b", "z")
    expect_equal(rownames(chrs), c("b", "z"))
})

test_that("featureData,fData,fData<-,MChromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    fd_exp <- annotatedDataFrameFrom(matrix(ncol = 2, nrow = 2), byrow = TRUE)
    expect_equal(featureData(chrs), fd_exp)

    ## Check error when assigning a featureData with different names
    fd <- data.frame(name = letters[1:2], idx = 1:2)
    rownames(fd) <- letters[1:2]
    expect_error(fData(chrs) <- pd)

    fd <- data.frame(name = letters[1:2], idx = 1:2)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          featureData = fd)
    expect_equal(featureData(chrs), AnnotatedDataFrame(fd))
    expect_equal(fData(chrs), fd)

    fd_2 <- cbind(fd, other = 1:2)
    fData(chrs) <- fd_2
    rownames(fd_2) <- as.character(1:2)
    expect_equal(featureData(chrs), AnnotatedDataFrame(fd_2))
    expect_equal(fData(chrs), fd_2)
    fd_3 <- cbind(fd_2, another = 3:4)
    featureData(chrs) <- fd_3
    expect_equal(featureData(chrs), AnnotatedDataFrame(fd_3))

    expect_equal(fvarLabels(chrs), colnames(fd_3))
})

test_that("isEmpty,MChromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    expect_true(!isEmpty(chrs))
    plot(chrs)

    chrs <- MChromatograms()
    expect_true(isEmpty(chrs))
    expect_warning(plot(chrs))

    ints <- rep(NA_real_, 105)
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- rep(NA_real_, 64)
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch1, ch2), nrow = 2)
    expect_true(isEmpty(chrs))
    expect_warning(plot(chrs))

    ## Only one row is empty.
    ints <- rep(NA_real_, 105)
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(64))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch1, ch2), nrow = 2)
    expect_true(!isEmpty(chrs))
    expect_warning(plot(chrs))

    ## 2x2 first row NA
    chrs <- MChromatograms(list(ch1, ch2, ch1, ch2), nrow = 2)
    expect_warning(plot(chrs))

    ## 2x2 first col NA
    chrs <- MChromatograms(list(ch1, ch1, ch2, ch2), nrow = 2)
    expect_true(!isEmpty(chrs))
    plot(chrs)
})

test_that(".mz_chromatograms, precursorMz etc,MChromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    ## Base function
    expect_error(MSnbase:::.mz_chromatograms(chrs, mz = "other"))
    chrs_f <- chrs
    fData(chrs_f) <- data.frame(precursorIsolationWindowTargetMZ = c(123, 456))
    res <- MSnbase:::.mz_chromatograms(chrs_f, mz = "precursorMz")
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs_f))
    expect_equal(res[, "mzmin"], res[, "mzmax"])
    expect_equal(res[, "mzmin"], c(123, 456))

    res <- MSnbase:::.mz_chromatograms(chrs_f)
    expect_true(all(is.na(res)))
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))

    ## method implementations:
    ## precursorMz
    expect_true(nrow(precursorMz(MChromatograms())) == 0)

    ## with precursor m/z data in the featureData data.frame
    chrs_f <- chrs
    fData(chrs_f) <- data.frame(precursorIsolationWindowTargetMZ = c(123, 456))
    res <- MSnbase:::.mz_chromatograms(chrs_f, "precursorMz")
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs_f))
    expect_equal(res[, "mzmin"], res[, "mzmax"])
    expect_equal(res[, "mzmin"], c(123, 456))
    expect_equal(precursorMz(chrs_f), res)

    ## Extracting precursor m/z data from the Chromatogram objects.
    res <- MSnbase:::.mz_chromatograms(chrs)
    expect_true(all(is.na(res)))
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))
    expect_equal(res, precursorMz(chrs))
    ## Define chrs.
    chrs_2 <- chrs
    chrs_2[1, 1]@precursorMz <- range(123)
    expect_error(MSnbase:::.mz_chromatograms(chrs_2, "precursorMz"),
                 "Chromatograms in row 1 have different precursorMz")
    chrs_2[1, 2]@precursorMz <- range(123)
    chrs_2[2, 1]@precursorMz <- range(456)
    chrs_2[2, 2]@precursorMz <- range(456)
    res <- MSnbase:::.mz_chromatograms(chrs_2, "precursorMz")
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))
    expect_equal(res[, "mzmin"], res[, "mzmax"])
    expect_equal(res[, "mzmin"], c(123, 456))
    expect_equal(res, precursorMz(chrs_2))

    ## productMz
    res <- MSnbase:::.mz_chromatograms(chrs, "productMz")
    expect_true(all(is.na(res)))
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))
    expect_equal(res, productMz(chrs))
    chrs_f <- chrs
    fData(chrs_f) <- data.frame(productIsolationWindowTargetMZ = c(3, 5))
    res <- MSnbase:::.mz_chromatograms(chrs_f, "productMz")
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))
    expect_equal(res[, "mzmin"], res[, "mzmax"])
    expect_equal(res[, "mzmin"], c(3, 5))
    expect_equal(res, productMz(chrs_f))

    chrs_2 <- chrs
    chrs_2[1, 1]@productMz <- range(5)
    expect_error(MSnbase:::.mz_chromatograms(chrs_2, "productMz"),
                 "Chromatograms in row 1 have different productMz")
    chrs_2[1, 2]@productMz <- range(5)
    res <- MSnbase:::.mz_chromatograms(chrs_2, "productMz")
    expect_equal(colnames(res), c("mzmin", "mzmax"))
    expect_equal(nrow(res), nrow(chrs))
    expect_equal(res[, "mzmin"], res[, "mzmax"])
    expect_equal(res[, "mzmin"], c(5, NA))
    expect_equal(res, productMz(chrs_2))

    ## polarity
    expect_true(all(polarity(chrs) == -1))
    fData(chrs)$polarity <- c(1, 1)
    expect_true(all(polarity(chrs) == 1))

    ## With a real object.
    on_disk <- microtofq_on_disk
    chrs <- chromatogram(on_disk, mz = c(123.4, 123.6), rt = c(35, 48))
    res <- mz(chrs)
    expect_true(nrow(res) == 1)
    expect_true(all(colnames(res) == c("mzmin", "mzmax")))
    expect_equal(unname(res[1, "mzmin"]), 123.4)
    expect_equal(unname(res[1, "mzmax"]), 123.6)
    expect_true(polarity(chrs) == unique(polarity(microtofq_on_disk)))
})

test_that(".bin_MChromatograms and bin,MChromatograms work", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- MChromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

    chrsb <- .bin_MChromatograms(chrs, binSize = 2)

    ## 1st row:
    expect_equal(rtime(chrsb[1, 1]), rtime(chrsb[1, 2]))
    expect_true(all(intensity(chrsb[1, 1])[rtime(chrsb[1, 1]) >
                                           max(rtime(chrs[1, 1]))] == 0))
    expect_true(max(rtime(chrsb[1, 2])) >= max(rtime(chrs[1, 2])))
    expect_true(max(rtime(chrsb[1, 1])) >= max(rtime(chrs[1, 1])))
    expect_equal(chrsb[1, 2], bin(chrs[1, 2], binSize = 2))
    ## 2nd row:
    expect_equal(rtime(chrsb[2, 1]), rtime(chrsb[2, 2]))
    expect_true(all(intensity(chrsb[2, 1])[rtime(chrsb[2, 1]) >
                                           max(rtime(chrs[2, 1]))] == 0))
    expect_true(max(rtime(chrsb[2, 2])) >= max(rtime(chrs[2, 2])))
    expect_true(max(rtime(chrsb[2, 1])) >= max(rtime(chrs[2, 1])))
    expect_equal(chrsb[2, 2], bin(chrs[2, 2], binSize = 2))
})

test_that("clean,MChromatograms works", {
    chr1 <- Chromatogram(rtime = 1:10,
                         intensity = c(0, 3, 5, 0, 0, 0, 8, 9, 0, 9))
    chr2 <- Chromatogram(rtime = 1:12,
                         intensity = c(0, 0, 0, 3, 5, 0, 0, 0, 8, 9, 0, 9))
    chr3 <- Chromatogram(rtime = 1:8,
                         intensity = c(3, 5, 0, 0, 8, 9, 0, 9))
    chr4 <- Chromatogram(rtime = 1:12,
                         intensity = c(NA, 5, 0, 0, 0, 8, 9, 0, 9, 0, 0, 0))
    chr5 <- Chromatogram(rtime = 1:5,
                         intensity = c(0, 0, 0, 0, 0))
    chr6 <- Chromatogram(rtime = 1:9,
                         intensity = c(0, 3, 5, 0, 0, 8, 9, 0, 9))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), nrow = 2,
                          phenoData = AnnotatedDataFrame(data.frame(id = 5:7)),
                          featureData = AnnotatedDataFrame(data.frame(mz = 3:4)))
    res <- clean(chrs)
    expect_equal(fData(res), fData(chrs))
    expect_equal(pData(res), pData(chrs))
    expect_equal(intensity(res[1, 1]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[2, 1]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[1, 2]), c(3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[2, 2]), c(NA, 5, 0, 0, 8, 9, 0, 9, 0))
    expect_equal(intensity(res[1, 3]), numeric())
    expect_equal(intensity(res[2, 3]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))

    res <- clean(chrs, na.rm = TRUE)
    expect_equal(intensity(res[1, 1]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[2, 1]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[1, 2]), c(3, 5, 0, 0, 8, 9, 0, 9))
    expect_equal(intensity(res[2, 2]), c(5, 0, 0, 8, 9, 0, 9, 0))
    expect_equal(intensity(res[1, 3]), numeric())
    expect_equal(intensity(res[2, 3]), c(0, 3, 5, 0, 0, 8, 9, 0, 9))

    res <- clean(chrs, all = TRUE)
    expect_equal(fData(res), fData(chrs))
    expect_equal(pData(res), pData(chrs))
    expect_equal(intensity(res[1, 1]), c(3, 5, 8, 9, 9))
    expect_equal(intensity(res[2, 1]), c(3, 5, 8, 9, 9))
    expect_equal(intensity(res[1, 2]), c(3, 5, 8, 9, 9))
    expect_equal(intensity(res[2, 2]), c(5, 8, 9, 9))
    expect_equal(intensity(res[1, 3]), numeric())
    expect_equal(intensity(res[2, 3]), c(3, 5, 8, 9, 9))
})

test_that("normalize,MChromatograms works", {
    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chr4 <- Chromatogram(rtime = 1:10,
                         intensity = c(NA, NA, 4, NA, NA, 9, NA, 10, 9, 1))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4), ncol = 2)
    res <- normalize(chrs)

    expect_true(ncol(res) == ncol(chrs))
    expect_true(nrow(res) == nrow(chrs))

    expect_equal(intensity(res[1, 2]) * max(intensity(chrs[1, 2]), na.rm = TRUE),
                 intensity(chrs[1, 2]))
})

test_that("filterIntensity,MChromatograms works", {
    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chr4 <- Chromatogram(rtime = 1:10,
                         intensity = c(NA, NA, 4, NA, NA, 9, NA, 10, 9, 1))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4), ncol = 2)
    res <- filterIntensity(chrs, intensity = 6)
    expect_true(ncol(res) == ncol(chrs))
    expect_true(nrow(res) == nrow(chrs))
    expect_equal(colnames(res), colnames(chrs))
    expect_equal(rownames(res), rownames(chrs))
    expect_true(all(intensity(res[1, 1]) >= 6))
    expect_true(all(intensity(res[1, 2]) >= 6))
    expect_true(all(intensity(res[2, 1]) >= 6))
    expect_true(all(intensity(res[2, 2]) >= 6))

    filt_fun <- function(x, prop = 0.2) {
        x@intensity >= max(x@intensity, na.rm = TRUE) * prop
    }

    res <- filterIntensity(chrs, intensity = filt_fun)
    expect_true(ncol(res) == ncol(chrs))
    expect_true(nrow(res) == nrow(chrs))
    expect_equal(colnames(res), colnames(chrs))
    expect_equal(rownames(res), rownames(chrs))
    expect_true(all(intensity(res[1, 1]) >=
                    max(intensity(chrs[1, 1]), na.rm = TRUE) * 0.2))
    expect_true(all(intensity(res[1, 2]) >=
                    max(intensity(chrs[1, 2]), na.rm = TRUE) * 0.2))
    expect_true(all(intensity(res[2, 1]) >=
                    max(intensity(chrs[2, 1]), na.rm = TRUE) * 0.2))
    expect_true(all(intensity(res[2, 2]) >=
                    max(intensity(chrs[2, 2]), na.rm = TRUE) * 0.2))

    res <- filterIntensity(chrs, intensity = filt_fun, prop = 0.5)
    expect_true(ncol(res) == ncol(chrs))
    expect_true(nrow(res) == nrow(chrs))
    expect_equal(colnames(res), colnames(chrs))
    expect_equal(rownames(res), rownames(chrs))
    expect_true(all(intensity(res[1, 1]) >=
                    max(intensity(chrs[1, 1]), na.rm = TRUE) * 0.5))
    expect_true(all(intensity(res[1, 2]) >=
                    max(intensity(chrs[1, 2]), na.rm = TRUE) * 0.5))
    expect_true(all(intensity(res[2, 1]) >=
                    max(intensity(chrs[2, 1]), na.rm = TRUE) * 0.5))
    expect_true(all(intensity(res[2, 2]) >=
                    max(intensity(chrs[2, 2]), na.rm = TRUE) * 0.5))
})

test_that(".compare_chromatograms works", {
    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    res <- .compare_chromatograms(list(chr1), list(chr2, chr3))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), 1)
    expect_true(res[1, 1] < 0.5)
    expect_true(res[1, 2] > 0.9)
    res2 <- .compare_chromatograms(list(chr1), list(chr2, chr3),
                                   ALIGNFUNARGS = list(tolerance = 0))
    expect_true(all(is.na(res2)))

    expect_equal(res, .compare_chromatograms(list(chr1), list(chr2, chr3)))

    res <- .compare_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2, chr3))
    expect_equal(res[1, 2], compareChromatograms(chr1, chr2))
    expect_equal(res[1, 3], compareChromatograms(chr1, chr3))
    expect_equal(res[2, 1], compareChromatograms(chr2, chr1))
    expect_equal(res[2, 3], compareChromatograms(chr2, chr3))
    expect_equal(res[3, 1], compareChromatograms(chr3, chr1))

    res <- .compare_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2, chr3))
    expect_equal(res[1, 2], compareChromatograms(chr1, chr2))
    expect_equal(res[1, 3], compareChromatograms(chr1, chr3))
    expect_equal(res[2, 1], res[1, 2])
    expect_equal(res[2, 3], compareChromatograms(chr2, chr3))
    expect_equal(res[3, 1], res[1, 3])

    res <- .compare_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2))

    chrs <- MChromatograms(list(chr1, chr2, chr3, chr1), ncol = 2)
    expect_error(.compare_chromatograms(chrs, list(chr1, chr2)),
                 "single column")
    expect_error(.compare_chromatograms(list(chr1, chr2), chrs),
                 "single column")
})

test_that("compareChromatograms,MChromatograms works", {
    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))

    res <- compareChromatograms(chrs)
    expect_true(nrow(res) == 3)
    expect_true(ncol(res) == 3)
    expect_true(res[1, 3] > 0.9)
    expect_true(res[1, 2] < 0.5)

    res_2 <- compareChromatograms(chrs, chrs)
    expect_equal(res_2, res)

    res_3 <- compareChromatograms(chrs, chrs[1:2, ])
    expect_true(nrow(res_3) == 3)
    expect_true(ncol(res_3) == 2)
    expect_equal(res_3, res_2[, 1:2])
})

test_that("transformIntensity works", {
    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chr4 <- Chromatogram(rtime = 1:10,
                         intensity = c(NA, NA, 4, NA, NA, 9, NA, 10, 9, 1))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4), ncol = 2)

    chrs2 <- transformIntensity(chrs)
    expect_equal(intensity(chrs[1, 1]), intensity(chrs2[1, 1]))
    expect_equal(intensity(chrs[1, 2]), intensity(chrs2[1, 2]))
    expect_equal(intensity(chrs[2, 1]), intensity(chrs2[2, 1]))
    expect_equal(intensity(chrs[2, 2]), intensity(chrs2[2, 2]))

    chrs2 <- transformIntensity(chrs, FUN = log2)
    expect_equal(log2(intensity(chrs[1, 1])), intensity(chrs2[1, 1]))
    expect_equal(log2(intensity(chrs[1, 2])), intensity(chrs2[1, 2]))
    expect_equal(log2(intensity(chrs[2, 1])), intensity(chrs2[2, 1]))
    expect_equal(log2(intensity(chrs[2, 2])), intensity(chrs2[2, 2]))
})
