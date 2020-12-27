test_that("imputation works", {
    data(naset)
    res1 <- impute(naset, method = "knn")
    res2 <- MsCoreUtils::impute_knn(exprs(naset))
    expect_identical(exprs(res1), res2)
})
