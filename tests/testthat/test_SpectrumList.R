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

test_that("SpectrumList construction works as expected", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    
    ## Errors.
    expect_error(new("SpectrumList", 4))
    expect_error(new("SpectrumList", list(4)))
    expect_error(new("SpectrumList", list(sp1, 4)))
    expect_error(SpectrumList(4))

    spl <- new("SpectrumList", list(sp1, sp2, sp3))
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)

    spl <- SpectrumList(sp1)
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    spl <- SpectrumList(sp1, sp2)
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    spl <- SpectrumList(list(sp1, sp2, sp3),
                        elementMetadata = DataFrame(id = c("a", "b", "c")))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)
    expect_equal(mcols(spl), DataFrame(id = c("a", "b", "c")))
    
    ## Concatenating.
    spl <- c(SpectrumList(sp1), SpectrumList(sp2))
    expect_true(is(spl, "SpectrumList"))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)

})

test_that(".make_naked_matrix_from_SpectrumList works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- SpectrumList(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))
    
    res <- .make_naked_matrix_from_SpectrumList(spl)
    expect_equal(ncol(res), 5)
    expect_equal(nrow(res), 3)
    
    spl <- SpectrumList(new("Spectrum2", mz = c(1, 2, 3), intensity = 1:3))
    res <- .make_naked_matrix_from_SpectrumList(spl)
    expect_equal(ncol(res), 3)
    expect_equal(nrow(res), 1)
})

test_that("show,SpectrumList works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- SpectrumList(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))

    .show_SpectrumList(spl)
    .show_SpectrumList(spl, print.classinfo = TRUE)
    show(spl)
})

