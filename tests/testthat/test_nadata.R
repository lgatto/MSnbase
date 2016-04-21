context("NA data")
library(pRolocdata)
data(dunkley2006)

test_that("makeNAdata", {
    expect_error(makeNaData(1))
    expect_error(makeNAdata(dunkley2006))
    tmp <- makeNaData(dunkley2006, nNA = 1, exclude = 1:688)

    k <- 689
    tmp <- makeNaData(dunkley2006, nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 1L)
    expect_identical(sum(is.na(exprs(tmp[-k, ]))), 0L)

    k <- logical(689)
    k[111] <- TRUE
    tmp <- makeNaData(dunkley2006, nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[!k, ]))), 1L)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 0L)

    k <- featureNames(dunkley2006)[222]
    tmp <- makeNaData(dunkley2006, nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 1L)
    expect_identical(sum(is.na(exprs(tmp[setdiff(featureNames(dunkley2006), k), ]))),
                     0L)

    tmp <- makeNaData(dunkley2006, nNA = 100)
    expect_identical(sum(is.na(exprs(tmp))), 100L)

    tmp <- makeNaData(dunkley2006, pNA = 0.1)
    expect_identical(sum(is.na(exprs(tmp))), 1103L)
})


test_that("makeNAdata2", {
    expect_error(makeNAdata2(1))
})

test_that("whichNA", {
    k <- 1
    tmp <- makeNaData(dunkley2006, nNA = 16, exclude = 2:689)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 16L)
    expect_identical(sum(is.na(exprs(tmp)[-k, ])), 0L)
})
