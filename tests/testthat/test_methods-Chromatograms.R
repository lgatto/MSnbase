test_that("[,Chromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    
    ## o Subset using indices
    expect_true(is(chrs[1, 1], "Chromatogram"))
    expect_equal(chrs[1, 2], ch2)
    ##   extract a row
    expect_equal(chrs[1, , drop = TRUE], list(`1` = ch, `2` = ch2))
    expect_equal(chrs[1, , drop = FALSE], Chromatograms(list(ch, ch2), nrow = 1))
    ##   Test the default
    expect_equal(chrs[1, ], Chromatograms(list(ch, ch2), nrow = 1))    
    ##   extract a column
    expect_equal(chrs[, 2, drop = TRUE], list(ch2, ch3))
    res <- chrs[, 2, drop = FALSE]
    res_exp <- Chromatograms(list(ch2, ch3), ncol = 1,
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
    res_exp <- Chromatograms(list(ch, ch2), nrow = 1)
    colnames(res_exp) <- c("a", "b")
    expect_equal(chrs[1, , drop = FALSE], res_exp)
    ##   Test the default
    expect_equal(chrs[1, ], res_exp)
    ##   extract a column
    expect_equal(chrs[, 2, drop = TRUE], list(ch2, ch3))
    res_exp <- Chromatograms(list(ch2, ch3), ncol = 1)
    colnames(res_exp) <- "b"
    expect_equal(chrs[, 2, drop = FALSE], res_exp)
    
    ## o Subset using logical
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    expect_true(is(chrs[c(TRUE, FALSE), c(TRUE, FALSE)], "Chromatogram"))
    expect_equal(chrs[c(TRUE, FALSE), c(FALSE, TRUE)], ch2)
    ##   extract a row
    expect_equal(chrs[c(TRUE, FALSE), , drop = TRUE], list(`1` = ch, `2` = ch2))
    expect_equal(chrs[c(TRUE, FALSE), , drop = FALSE],
                 Chromatograms(list(ch, ch2), nrow = 1))
    expect_equal(chrs[c(TRUE, FALSE), ],
                 Chromatograms(list(ch, ch2), nrow = 1))
    ##   extract a column
    expect_equal(chrs[, c(FALSE, TRUE), drop = TRUE], list(ch2, ch3))
    res <- chrs[, c(FALSE, TRUE), drop = FALSE]
    rownames(pData(res)) <- rownames(pData(res))
    expect_equal(res, Chromatograms(list(ch2, ch3), ncol = 1,
                                    dimnames = list(NULL, "2")))
    ##   Repeat with colnames
    colnames(chrs) <- c("a", "b")
    expect_equal(chrs[c(TRUE, FALSE), , drop = TRUE], list(a = ch, b = ch2))
    res_exp <- Chromatograms(list(ch, ch2), nrow = 1)
    colnames(res_exp) <- c("a", "b")
    expect_equal(chrs[c(TRUE, FALSE), , drop = FALSE], res_exp)
    expect_equal(chrs[c(TRUE, FALSE), ], res_exp)
    ## extract a column
    expect_equal(chrs[, c(FALSE, TRUE), drop = TRUE], list(ch2, ch3))
    res_exp <- Chromatograms(list(ch2, ch3), ncol = 1)
    colnames(res_exp) <- "b"
    expect_equal(chrs[, c(FALSE, TRUE)], res_exp)
    expect_equal(chrs[, c(FALSE, TRUE)], res_exp)

    ## Subset using names
    expect_equal(chrs[, "a", drop = TRUE], list(ch, ch1))
    res_exp <- Chromatograms(list(ch, ch1), ncol = 1)
    colnames(res_exp) <- "a"
    expect_equal(chrs[, "a", drop = FALSE], res_exp)
    expect_equal(chrs[, "a"], res_exp)

    ## Check phenoData while subsetting.
    pd <- data.frame(name = letters[1:2], idx = 1:2)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    res <- chrs[, 2]
    pd_exp <- droplevels(pd[2, ])
    expect_equal(pData(res), pd_exp)
    rownames(pd) <- c("g", "h")
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    res <- chrs[, 2]
    pd_exp <- droplevels(pd[2, ])
    expect_equal(pData(res), pd_exp)
})

test_that("[<-,Chromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
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
    expect_equal(chrs[, 2, drop = TRUE], list(ch2, ch3))
    chrs[2, 1] <- list(ch4)
    expect_equal(chrs[2, 1], ch4)
    
    chrs[, "a"] <- list(ch2, ch1)
    expect_equal(chrs[, 1, drop = TRUE], list(ch2, ch1))
    expect_error(chrs[, 1] <- list(ch, ch2, ch3))
    
    chrs[, c(TRUE, FALSE)] <- list(ch4, ch4)
    expect_equal(chrs[, 1, drop = TRUE], list(ch4, ch4))
    
})

test_that("plot,Chromatograms works", {
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
    plot(chrs)
    plot(chrs[1, , drop = FALSE])
    plot(chrs[1, 1, drop = FALSE])
    plot(chrs[1, ])
    plot(chrs[1, 1])
    plot(chrs[, 2])
})

test_that("colnames<-, sampleNames, sampleNames<-,Chromatograms works", {
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)

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

test_that("phenoData,pData,pData<-,Chromatograms works", {
    ## Check if we can access the phenoData.
    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    pd_exp <- annotatedDataFrameFrom(matrix(ncol = 2, nrow = 2), byrow = FALSE)
    rownames(pData(pd_exp)) <- NULL
    pd_exp <- as(pd_exp, "NAnnotatedDataFrame")
    expect_equal(phenoData(chrs), pd_exp)
    
    pd <- data.frame(name = letters[1:2], idx = 1:2)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                          phenoData = AnnotatedDataFrame(pd))
    chrs_2 <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2,
                            phenoData = pd)
    expect_equal(chrs, chrs_2)
    expect_equal(phenoData(chrs), as(AnnotatedDataFrame(pd), "NAnnotatedDataFrame"))

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
    expect_equal(chrs$name, factor(letters[1:2]))
    expect_equal(chrs$idx, 1:2)

    chrs$idx <- c(2, 1)
    expect_equal(chrs$idx, c(2, 1))

    chrs$new_variable <- c("it", "works")
    expect_equal(chrs$new_variable, c("it", "works"))

    expect_error(chrs$new_variable <- 1:4)
    chrs$new_variable <- 1
})

