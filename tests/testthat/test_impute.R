context("imputation method")

test_that("all imputation methods", {

    data(naset)
  
    m <- imputeMethods()
    m <- m[m != "mixed"]
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
