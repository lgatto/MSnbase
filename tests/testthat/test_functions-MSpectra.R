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

test_that("MSpectra construction works as expected", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)

    ## Errors.
    expect_error(new("MSpectra", 4))
    expect_error(new("MSpectra", list(4)))
    expect_error(new("MSpectra", list(sp1, 4)))
    expect_error(MSpectra(4))

    spl <- new("MSpectra", list(sp1, sp2, sp3))
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)

    spl <- MSpectra(sp1)
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    spl <- MSpectra(sp1, sp2)
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    spl <- MSpectra(list(sp1, sp2, sp3),
                   elementMetadata = DataFrame(id = c("a", "b", "c")))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    expect_equal(spl[[3]], sp3)
    expect_equal(mcols(spl), DataFrame(id = c("a", "b", "c"),
                                       row.names = 1:3))

    ## Concatenating.
    spl <- c(MSpectra(sp1), MSpectra(sp2))
    expect_true(is(spl, "MSpectra"))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)

})

test_that(".make_naked_matrix_from_MSpectra works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- MSpectra(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))

    res <- .make_naked_matrix_from_MSpectra(spl)
    expect_equal(ncol(res), 5)
    expect_equal(nrow(res), 3)

    spl <- MSpectra(new("Spectrum2", mz = c(1, 2, 3), intensity = 1:3))
    res <- .make_naked_matrix_from_MSpectra(spl)
    expect_equal(ncol(res), 3)
    expect_equal(nrow(res), 1)
})

test_that("show,MSpectra works", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
               precursorMz = 2, rt = 1.232446)
    sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
               rt = 1.232445)
    spl <- MSpectra(sp1, sp2, sp3, elementMetadata = DataFrame(id = 1:3))

    .show_MSpectra(spl)
    .show_MSpectra(spl, print.classinfo = TRUE)
    show(spl)
})

test_that("extractSpectraData works", {
    fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
                      package = "msdata")
    data <- filterRt(readMSData(fl, mode = "onDisk"), rt = c(1, 6))
    sps <- spectra(data)

    res <- extractSpectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_true(all(c("mz", "intensity") %in% colnames(res)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    expect_error(extractSpectraData(1:10), "should be either a 'list'")

    res <- extractSpectraData(data)
    expect_true(is(res, "DataFrame"))
    expect_true(all(c("mz", "intensity") %in% colnames(res)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    spctra <- MSpectra(sps)
    mcols(spctra)$new_col <- "a"
    res <- extractSpectraData(sps)
    expect_true(all(res$new_col == "a"))
})
