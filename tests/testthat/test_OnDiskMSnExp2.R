context("OnDiskMSnExp class")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
x1 <- readMSData(f, verbose = FALSE)
x2 <- MSnbase:::readMSData2(f, msLevel = 2)

test_that("compare MS2 on disk and in memoery", {
    expect_identical(length(x1), length(x2))
    expect_false(identical(featureNames(x1), featureNames(x2)))
    featureNames(x2) <- featureNames(x1)

})
