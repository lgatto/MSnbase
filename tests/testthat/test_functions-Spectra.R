test_that(".short_spectrum_info works", {
    sp1 <- new("Spectrum1", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 3)
    expect_equal(.short_spectrum_info(sp1), c(msLevel = sp1@msLevel,
                                              rtime = NA,
                                              peaksCount = peaksCount(sp1)))
    expect_equal(.short_spectrum_info(sp2), c(msLevel = sp2@msLevel,
                                              rtime = sp2@rt,
                                              peaksCount = peaksCount(sp2)))
})

test_that("Spectra construction works as expected", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)

    ## Errors.
    expect_error(new("Spectra", 4))
    expect_error(new("Spectra", list(4)))
    expect_error(new("Spectra", list(sp1, 4)))
    expect_error(Spectra(4))

    spl <- new("Spectra", list(sp1, sp2, sp3))
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)

    spl <- Spectra(sp1)
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    spl <- Spectra(sp1, sp2)
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    spl <- Spectra(list(sp1, sp2, sp3),
                   elementMetadata = DataFrame(id = c("a", "b", "c")))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)
    expect_equal(mcols(spl), DataFrame(id = c("a", "b", "c"),
                                       row.names = 1:3))

    ## Concatenating.
    spl <- c(Spectra(sp1), Spectra(sp2))
    expect_true(is(spl, "Spectra"))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)

})

test_that(".make_naked_matrix_from_Spectra works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- Spectra(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))

    res <- .make_naked_matrix_from_Spectra(spl)
    expect_equal(ncol(res), 5)
    expect_equal(nrow(res), 3)

    spl <- Spectra(new("Spectrum2", mz = c(1, 2, 3), intensity = 1:3))
    res <- .make_naked_matrix_from_Spectra(spl)
    expect_equal(ncol(res), 3)
    expect_equal(nrow(res), 1)
})

test_that("show,Spectra works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- Spectra(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))

    .show_Spectra(spl)
    .show_Spectra(spl, print.classinfo = TRUE)
    show(spl)
})

test_that(".combine_spectra works", {
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)
    res <- .combine_spectra(sps)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], sort(unlist(spd$mz)))

    res <- .combine_spectra(sps, FUN = combinePeaks, tolerance = 0.1)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(mean(c(12, 12.1)), mean(c(14, 14.1, 14.15)),
                                mean(c(34, 34.1)), 45, mean(c(56, 56.1))))
    expect_equal(res$intensity[[1]], c(mean(c(10, 12)), mean(c(20, 11, 22)),
                                       mean(c(21, 32)), 30, mean(c(40, 31))))
    res <- .combine_spectra(sps, FUN = combinePeaks, tolerance = 0.1,
                            mzFun = max, intensityFun = median)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(max(c(12, 12.1)), max(c(14, 14.1, 14.15)),
                                max(c(34, 34.1)), 45, max(c(56, 56.1))))
    expect_equal(res$intensity[[1]], c(median(c(10, 12)), median(c(20, 11, 22)),
                                       median(c(21, 32)), 30, median(c(40, 31))))
})
