context("imputation method")

test_that("all imputation methods", {
    data(naset, package = "MSnbase")
    m <- imputeMethods()
    m <- m[m != "mixed"]
    m <- m[m != "nbavg"] ## see next test
    
    for (.m in m) {
        xx <- impute(naset, method = .m)
        expect_true(validObject(xx))
        expect_false(any(is.na(exprs(xx))))
    }

    expect_error(impute(naset, method = "mixed",
                        randna = fData(naset)$randna,
                        mnar = "min"),
                 regexp = "mar")

    expect_error(impute(naset, method = "mixed",
                        randna = fData(naset)$randna,
                        mar = "knn"),
                 regexp = "mnar")

    expect_error(impute(naset, method = "mixed",
                        mnar = "min",
                        mar = "knn"),
                 regexp = "randna")

    expect_error(impute(naset, method = "mixed",
                        randna = TRUE,
                        mnar = "min",
                        mar = "knn"),
                 regexp = "randna")

    mx <- impute(naset, method = "mixed",
                 randna = fData(naset)$randna,
                 mnar = "min",
                 mar = "knn")
    
    expect_true(validObject(mx))
    expect_false(any(is.na(exprs(mx))))
})


test_that("nbavg methods", {
    
    m <- matrix(1:25, 5)
    ## default min value
    m[1, 2] <- 0.1 
    ## imputes as min value (or use-defined k)
    m[1, 1] <- m[5, 5] <- NA 
    m[2, 1:2] <- NA ## [2, 1] will be min
                    ## [2, 2] will be avg 6.05
    ## remaing NA    
    m[3, 3:4] <- NA
    ## average imputation
    m[5, 2] <- NA ## will be 10
    m[4, 3] <- NA ## will be 14    
    pd <- fd <- data.frame(A = 1:5)
    rownames(m) <- colnames(m) <-
        rownames(pd) <- rownames(fd) <-
            LETTERS[1:5]
    x <- MSnSet(m, fd, pd)

    xx <- impute(x, "nbavg")
    expect_true(exprs(xx[1, 2]) == 0.1)
    expect_true(exprs(xx[1, 1]) == 0.1)
    expect_true(exprs(xx[2, 1]) == 0.1)
    expect_true(exprs(xx[2, 2]) == 6.05)
    expect_true(all(is.na(exprs(xx[3, 3:4]))))
    expect_true(exprs(xx[5, 2]) == 10)
    expect_true(exprs(xx[4, 3]) == 14)

    xx <- impute(x, "nbavg", k = 0)
    expect_true(exprs(xx[1, 2]) == 0.1)
    expect_true(exprs(xx[1, 1]) == 0)
    expect_true(exprs(xx[2, 1]) == 0)
    expect_true(exprs(xx[2, 2]) == 6)
    expect_true(all(is.na(exprs(xx[3, 3:4]))))
    expect_true(exprs(xx[5, 2]) == 10)
    expect_true(exprs(xx[4, 3]) == 14)
})
