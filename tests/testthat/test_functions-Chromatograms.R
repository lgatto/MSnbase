test_that("Chromatograms works", {
    chs <- new("Chromatograms")
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    chs <- Chromatograms()
    expect_equal(nrow(chs), 0)
    expect_equal(ncol(chs), 0)
    ## Errors:
    chs@.Data <- matrix(1:4)
    expect_error(validObject(chs))
    expect_error(Chromatograms(1:4))
    ch <- new("Chromatogram")
    ch_list <- list(ch, ch, ch, ch, ch, ch, ch, ch)
    chs <- Chromatograms(ch_list, nrow = 2)
    expect_equal(length(chs[1, 1]), 0)
    chs[2, 1] <- list(Chromatogram(1:10, 1:10))
    expect_equal(length(chs[2, 1]), 10)
    expect_equal(chs[, 1], c(ch, Chromatogram(1:10, 1:10)))
})

test_that(".plotChromatogramList works", {
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
    MSnbase:::.plotChromatogramList(chrs[1, ])
    MSnbase:::.plotChromatogramList(chrs[1, ], main = "some title",
                                    col = c("red", "blue"))
    expect_error(MSnbase:::.plotChromatogramList(1:10))
})
