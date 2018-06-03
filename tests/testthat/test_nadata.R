context("NA data")
library("pRolocdata")
data(dunkley2006)

test_that("makeNAdata", {
    expect_error(makeNaData(1))
    expect_error(makeNAdata(dunkley2006))
    tmp <- makeNaData(dunkley2006, nNA = 1, exclude = 1:688)

    k <- 1
    tmp <- makeNaData(dunkley2006[1:2, ], nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 0L)
    expect_identical(sum(is.na(exprs(tmp[-k, ]))), 1L)

    k <- logical(2)
    k[1] <- TRUE
    tmp <- makeNaData(dunkley2006[1:2, ], nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[!k, ]))), 1L)
    expect_identical(sum(is.na(exprs(tmp[k, ]))), 0L)

    k <- featureNames(dunkley2006)[1]
    tmp <- makeNaData(dunkley2006[1:2, ], nNA = 1, exclude = k)
    expect_identical(sum(is.na(exprs(tmp[1, ]))), 0L)
    expect_identical(sum(is.na(exprs(tmp[2, ]))), 1L)

    tmp <- makeNaData(dunkley2006, nNA = 100)
    expect_identical(sum(is.na(exprs(tmp))), 100L)

    tmp <- makeNaData(dunkley2006, pNA = 0.1)
    expect_identical(sum(is.na(exprs(tmp))), 1103L)
})


test_that("makeNAdata2", {
    expect_error(makeNaData2(1))
    expect_error(makeNaData2(dunkley2006))
    expect_error(makeNaData2(dunkley2006, nRows = 1:10, nNA = 1:2))


    tmp <- makeNaData2(dunkley2006[1:56],
                       nRows = 1:10, nNA = 1:10,
                       exclude = 1)
    expect_true(!anyNA(tmp[1, ]))
    expect_true(all(apply(exprs(tmp[-1, ]), 1, function(x) any(is.na(x)))))
    nna <- table(apply(exprs(tmp[-1, ]), 1, function(x) sum(is.na(x))))
    expect_identical(as.vector(nna), 1:10)

    k <- logical(56)
    k[1] <- TRUE
    tmp <- makeNaData2(dunkley2006[1:56],
                       nRows = 1:10, nNA = 1:10,
                       exclude = k)
    expect_true(!anyNA(tmp[1, ]))
    expect_true(all(apply(exprs(tmp[-1, ]), 1, function(x) any(is.na(x)))))
    nna <- table(apply(exprs(tmp[-1, ]), 1, function(x) sum(is.na(x))))
    expect_identical(as.vector(nna), 1:10)

    k <- featureNames(dunkley2006)[1]
    tmp <- makeNaData2(dunkley2006[1:56],
                       nRows = 1:10, nNA = 1:10,
                       exclude = k)
    expect_true(!anyNA(tmp[1, ]))
    expect_true(all(apply(exprs(tmp[-1, ]), 1, function(x) any(is.na(x)))))
    nna <- table(apply(exprs(tmp[-1, ]), 1, function(x) sum(is.na(x))))
    expect_identical(as.vector(nna), 1:10)

    tmp <- makeNaData2(dunkley2006, nRows = 1:10, nNA = 1:10)
    expect_identical(sum(is.na(exprs(tmp))), sum(1:10 * 1:10))
})

test_that("whichNA", {
    k <- 1
    tmp <- makeNaData(dunkley2006, nNA = 16, exclude = 2:689)
    wna <- whichNA(tmp)
    expect_identical(matrix(c(rep(1L, 16), 1:16), ncol = 2), wna)
})
