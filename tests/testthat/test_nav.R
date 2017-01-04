context("Navigating MSnExp data")

library("msdata")
f <- proteomics(full.names = TRUE)[1]


test_that("nextMS and prevMS", {
    x <- readMSData2(f)
    sp <- "X009.1"
    expect_identical(acquisitionNum(x[[sp]]), 32926L)
    expect_identical(acquisitionNum(MSnbase:::nextMS(sp, x, 3L)), 32929L)
    expect_identical(acquisitionNum(MSnbase:::prevMS(sp, x, 3L)), 32923L)
    expect_identical(acquisitionNum(MSnbase:::nextMS(sp, x, 2L)), 32927L)
    expect_identical(acquisitionNum(MSnbase:::prevMS(sp, x, 2L)), 32925L)
    expect_identical(acquisitionNum(MSnbase:::nextMS(sp, x, 1L)), 32941L)
    expect_identical(acquisitionNum(MSnbase:::prevMS(sp, x, 1L)), 32918L)
    expect_null(MSnbase:::prevMS("X001.1", x, 1L))
    expect_null(MSnbase:::prevMS("X001.1", x, 2L))
    expect_null(MSnbase:::prevMS("X001.1", x, 3L))
    expect_identical(MSnbase:::nextMS(sp, x, 3L), MSnbase:::nextMS(sp, x))
})
