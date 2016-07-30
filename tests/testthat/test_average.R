library("pRolocdata")

test_that("averageMSnSet", {
    m1 <- matrix(1, ncol = 4, nrow = 4)
    m2 <- matrix(2, ncol = 4, nrow = 4)
    m3 <- matrix(3, ncol = 4, nrow = 4)
    m2[1, 1] <- m1[1, 1] <- NA
    m3[4, 1] <- 10
    fd1 <- fd2 <- fd3 <- data.frame(X = 1:4)
    rownames(fd1) <- rownames(m1) <- letters[1:4]
    rownames(fd2) <- rownames(m2) <- letters[2:5]
    rownames(fd3) <- rownames(m3) <- letters[3:6]
    colnames(m1) <- colnames(m2)  <- LETTERS[1:4]
    colnames(m3) <- LETTERS[4:1]
    pd <- data.frame(s = 1:4, row.names = colnames(m1))
    pd3 <- data.frame(s = 1:4, row.names = colnames(m3))
    e1 <- MSnSet(m1, fd1, pd)
    e2 <- MSnSet(m2, fd2, pd)
    e3 <- MSnSet(m3, fd3, pd3)
    expect_true(validObject(e1))
    expect_true(validObject(e2))
    expect_true(validObject(e3))

    le <- MSnSetList(list(e1, e2, e3))
    avg <- averageMSnSet(le)
    expect_equal(dim(avg), c(6, 4))
    expect_equal(sampleNames(e1), sampleNames(avg))
    expect_equal(featureNames(avg), letters[1:6])

    res <- matrix(c(NaN, rep(1.0, 3),
                    1, rep((1+2)/2, 3),
                    rep((1+2+3)/3, 8),
                    rep((2+3)/2, 4),
                    rep(3, 3), 10),
                  byrow = TRUE,
                  ncol = 4, nrow = 6)
    colnames(res) <- LETTERS[1:4]
    rownames(res) <- letters[1:6]

    expect_equal(res, exprs(avg))
    expect_equal(dim(avg), dim(fData(avg)$nNA))
    expect_equal(dim(avg), dim(fData(avg)$disp))

    fd <- t(data.frame(a = c(3, 2, 2, 2),
                       b = c(2, 1, 1, 1),
                       c = rep(0, 4),
                       d = rep(0, 4),
                       e = rep(1, 4),
                       f = rep(2, 4),
                       row.names = LETTERS[1:4]))
    expect_equal(fData(avg)$nNA, fd)

    ## using sd would be
    nadisp <- t(data.frame(a = rep(TRUE, 4),
                           b = c(TRUE, rep(FALSE, 3)),
                           c = rep(FALSE, 4),
                           d = rep(FALSE, 4),
                           e = rep(FALSE, 4),
                           f = rep(TRUE, 4),
                           row.names = LETTERS[1:4]))
    ## using mad
    nadisp <- matrix(FALSE, ncol = 4, nrow = 6)
    rownames(nadisp) <- letters[1:6]
    colnames(nadisp) <- LETTERS[1:4]
    nadisp[1, 1] <- TRUE

    expect_equal(is.na(fData(avg)$disp),
                 nadisp)
    sampleNames(e3) <- LETTERS[10:13]
    expect_error(averageMSnSet(MSnSetList(list(e1, e3))))
})

test_that("averageMSnSet list of length 1", {
    data(dunkley2006, package = "pRolocdata")
    expect_equal(dunkley2006,
                 averageMSnSet(MSnSetList(list(dunkley2006))))
})
