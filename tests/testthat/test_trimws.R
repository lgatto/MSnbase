context("trimming white spaces")

## make a test data
library("pRolocdata")
data(dunkley2006)
fData(dunkley2006)$test0 <-
                     fData(dunkley2006)$test <-
                                          as.character(1:nrow(dunkley2006))
fData(dunkley2006)$test[1:10] <- paste0(" ", fData(dunkley2006)$test[1:10])
fData(dunkley2006)$test[11:20] <- paste0(fData(dunkley2006)$test[11:20], " ")
fData(dunkley2006)$test[21:30] <- paste0("  ", fData(dunkley2006)$test[21:30], " ")


test_that("trimws on MSnSet", {
    d2 <- trimws(dunkley2006)
    expect_identical(fData(d2)$test, fData(d2)$test0)

    d2 <- trimws(dunkley2006, which = "left")
    expect_false(identical(fData(d2)$test, fData(d2)$test0))
    expect_false(identical(fData(d2)$test[-(1:10)], fData(d2)$test0[-(1:10)]))
    expect_identical(fData(d2)$test[1:10], fData(d2)$test0[1:10])

    d2 <- trimws(dunkley2006, which = "right")
    expect_false(identical(fData(d2)$test, fData(d2)$test0))
    expect_false(identical(fData(d2)$test[-(11:20)], fData(d2)$test0[-(11:20)]))
    expect_identical(fData(d2)$test[11:20], fData(d2)$test0[11:20])
})


test_that("trimws data.frame vs MSnSet ", {
    d2 <- trimws(dunkley2006)
    dfr2 <- trimws(fData(dunkley2006))
    expect_identical(fData(d2), dfr2)
})

test_that("check that base:::trimws still works", {
    expect_equal(trimws("  xx  "), "xx")
    expect_equal(trimws("  xx  ", "left"), "xx  ")
    expect_equal(trimws("  xx  ", "right"), "  xx")
})
