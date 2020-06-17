context("MSpectra-methods")

sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
           precursorMz = 2, rt = 1.232446)
sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
           rt = 1.232445)
spl_ <- MSpectra(sp1, sp2, sp3)

test_that("mz, intensity, rtime work", {
    spl <- spl_
    expect_true(length(mz(spl)) == 3)
    expect_equal(mz(spl[[2]]), mz(spl)[[2]])
    expect_equal(mz(spl)[[3]], mz(sp3))

    expect_true(length(intensity(spl)) == 3)
    expect_equal(intensity(spl[[2]]), intensity(spl)[[2]])
    expect_equal(intensity(spl)[[3]], intensity(sp3))

    expect_true(length(rtime(spl)) == 3)
    expect_equal(rtime(spl[[2]]), rtime(spl)[[2]])
    expect_equal(rtime(spl)[[3]], rtime(sp3))

    spl <- c(spl, MSpectra(new("Spectrum2")))
    expect_true(lengths(mz(spl))[4] == 0)
    expect_true(lengths(intensity(spl))[4] == 0)
    expect_equal(rtime(spl)[[4]], NA_real_)

    ## Put names on it.
    names(spl) <- c("a", "b", "c", "d")
    expect_equal(names(rtime(spl)), c("a", "b", "c", "d"))
    expect_equal(names(mz(spl)), c("a", "b", "c", "d"))
    expect_equal(names(intensity(spl)), c("a", "b", "c", "d"))

    ## Empty spectra
    spl <- MSpectra(new("Spectrum1"), new("Spectrum2"))
    expect_equal(rtime(spl), c(`1` = NA_real_, `2` = NA_real_))
    expect_true(length(mz(spl)) == 2)
    expect_true(all(lengths(mz(spl)) == 0))
    expect_true(length(intensity(spl)) == 2)
    expect_true(all(lengths(intensity(spl)) == 0))
})

test_that("precursor* work", {
    sp1 <- new("Spectrum2", precursorMz = 123.3, precursorCharge = 1L,
               precursorIntensity = 1234.4)
    sp2 <- new("Spectrum2", precursorMz = NA_real_, precursorCharge = integer(),
               precursorIntensity = NA_real_, precScanNum = 34L)
    spl <- MSpectra(sp1, sp2, sp3)
    expect_equal(precursorMz(spl), c(`1` = 123.3, `2` = NA, `3` = NA))
    expect_equal(precursorCharge(spl), c(`1` = 1L, `2` = NA, `3` = NA))
    expect_true(is.integer(precursorCharge(spl)))
    expect_equal(precursorIntensity(spl), c(`1` = 1234.4,`2` =  NA, `3` = NA))
    expect_equal(precScanNum(spl), c(`1` = NA_integer_, `2` = 34L, `3` = NA_integer_))
    expect_true(is.integer(precScanNum(spl)))

    expect_equal(precursorMz(spl_), c(`1` = NA_real_, `2` = 2, `3` = NA_real_))
    expect_equal(unname(precursorCharge(spl_)), rep(NA_integer_, length(spl_)))
    expect_equal(unname(precursorIntensity(spl_)), rep(NA_real_, length(spl_)))
    expect_true(is.integer(precScanNum(spl_)))
})

test_that("acquisitionNum and scanIndex work", {
    sp1 <- new("Spectrum2", acquisitionNum = 2L, scanIndex = 1L)
    sp2 <- new("Spectrum2", acquisitionNum = 4L)
    spl <- MSpectra(sp1, sp2)
    expect_identical(acquisitionNum(spl), c(`1` = 2L, `2` = 4L))
    expect_identical(scanIndex(spl), c(`1` = 1L, `2` = NA_integer_))

    expect_equal(unname(acquisitionNum(spl_)), rep(NA_integer_, length(spl_)))
    expect_equal(unname(scanIndex(spl_)), rep(NA_integer_, length(spl_)))
    expect_true(is.integer(acquisitionNum(spl_)))
    expect_true(is.integer(scanIndex(spl_)))
})

test_that("peaksCount, msLevel, tic and ionCount work", {
    sp1 <- new("Spectrum2", msLevel = 3L, tic = 5)
    sp2 <- new("Spectrum2")
    spl <- MSpectra(sp1, sp2)

    expect_true(is.integer(peaksCount(spl)))
    expect_equal(peaksCount(spl), c(`1` = 0, `2` = 0))
    expect_true(is.integer(msLevel(spl)))
    expect_equal(msLevel(spl), c(`1` = 3, `2` = 2))
    expect_true(is.numeric(tic(spl)))
    expect_equal(tic(spl), c(`1` = 5, `2` = 0))
    expect_true(is.numeric(ionCount(spl)))
    expect_equal(ionCount(spl), c(`1` = 0, `2` = 0))

    expect_equal(peaksCount(spl_), c(`1` = 3, `2` = 4, `3` = 5))
    expect_equal(msLevel(spl_), c(`1` = 2, `2` = 2, `3` = 1))
    expect_equal(tic(spl_), c(`1` = 11, `2` = 15, `3` = 40))
    expect_equal(ionCount(spl_), unlist(lapply(intensity(spl_), sum)))
})

test_that("collisionEnergy works", {
    sp1 <- new("Spectrum2")
    sp2 <- new("Spectrum2", collisionEnergy = 23.3)
    spl <- MSpectra(sp1, sp2)

    expect_true(is.numeric(collisionEnergy(spl)))
    expect_equal(collisionEnergy(spl), c(`1` = NA, `2` = 23.3))

    expect_true(is.numeric(collisionEnergy(spl_)))
    expect_equal(collisionEnergy(spl_), c(`1` = NA_real_, `2` = NA_real_, `3` = NA_real_))
})

test_that("fromFile and polarity work", {
    sp1 <- new("Spectrum2", polarity = 1L, fromFile = 5L)
    sp2 <- new("Spectrum2", fromFile = 3L)
    spl <- MSpectra(sp1, sp2)

    expect_true(is.integer(fromFile(spl)))
    expect_equal(fromFile(spl), c(`1` = 5, `2` = 3))
    expect_true(is.integer(polarity(spl)))
    expect_equal(polarity(spl), c(`1` = 1, `2` = NA))

    expect_true(is.integer(fromFile(spl_)))
    expect_true(all(is.na(fromFile(spl_))))
    expect_true(is.integer(polarity(spl_)))
    expect_equal(unname(polarity(spl_)), rep(NA_integer_, 3))
})

test_that("smoothed, isEmpty, centroided and isCentroided work", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(4, 2, 4, 5))
    sp2 <- new("Spectrum2", centroided = TRUE, smoothed = TRUE)
    spl <- MSpectra(sp1, sp2)

    expect_true(is.logical(smoothed(spl)))
    expect_equal(smoothed(spl), c(`1` = NA, `2` = TRUE))
    expect_true(is.logical(isEmpty(spl)))
    expect_equal(isEmpty(spl), c(`1` = FALSE, `2` = TRUE))
    expect_true(is.logical(centroided(spl)))
    expect_equal(centroided(spl), c(`1` = NA, `2` = TRUE))
    expect_true(is.logical(isCentroided(spl)))
    expect_equal(isCentroided(spl), c(`1` = NA, `2` = NA))

    expect_true(is.logical(smoothed(spl_)))
    expect_equal(unname(smoothed(spl_)), rep(NA, length(spl_)))
    expect_true(is.logical(isEmpty(spl_)))
    expect_true(all(!isEmpty(spl_)))
    expect_true(is.logical(centroided(spl_)))
    expect_true(all(is.na(centroided(spl_))))
    expect_true(is.logical(isCentroided(spl_)))
    expect_equal(isCentroided(spl_), c(`1` = NA, `2` = NA, `3` = NA))
})

test_that("writeMgfData,MSpectra works", {
    tmpf <- tempfile()

    writeMgfData(spl_, tmpf)
    res <- readLines(tmpf)
    ## No additional fields.
    expect_equal(res[7], "1 4")
    expect_equal(res[17], "1 5")
    expect_equal(res[27], "1 6")

    expect_error(writeMgfData(spl_, tmpf))

    spl <- spl_
    mcols(spl) <- DataFrame(index = 1:3, some_id = c("sp_1", "sp_2", "sp_3"))
    file.remove(tmpf)

    writeMgfData(spl, tmpf)
    res_2 <- readLines(tmpf)
    expect_equal(res_2[7], "INDEX=1")
    expect_equal(res_2[8], "SOME_ID=sp_1")
    expect_equal(res_2[19], "INDEX=2")
    expect_equal(res_2[20], "SOME_ID=sp_2")
    expect_equal(res_2[31], "INDEX=3")
    expect_equal(res_2[32], "SOME_ID=sp_3")
    expect_equal(res[-1], res_2[-c(1, 7, 8, 19, 20, 31, 32)])
})

test_that("clean,MSpectra works", {
    expect_equal(spl_, clean(spl_))
    spl <- c(spl_, MSpectra(new("Spectrum2", intensity = c(3, 0, 0, 3, 0, 5),
                               mz = c(1, 2, 3, 4, 5, 6))))
    res <- clean(spl, all = TRUE)
    expect_true(is(res, "MSpectra"))
    expect_equal(res[[4]], new("Spectrum2", intensity = c(3, 3, 5),
                               mz = c(1, 4, 6)))

    msnexp_clnd <- clean(tmt_od_sub)
    res_spctra <- clean(MSpectra(spectra(tmt_od_sub)))
    expect_equal(spectra(msnexp_clnd), res_spctra@listData)
})

test_that("removePeaks,MSpectra works", {
    spl <- MSpectra(new("Spectrum2", mz = 1:5, intensity = c(4, 34, 5, 199, 3),
                       centroided = TRUE),
                   new("Spectrum2", mz = 1:4, intensity = c(123, 43, 3, 9),
                       centroided = TRUE),
                   new("Spectrum1", mz = 1:4, intensity = c(123, 43, 2, 9),
                       centroided = TRUE),
                   elementMetadata = DataFrame(id = 1:3))
    res <- removePeaks(spl)
    expect_true(is(res, "MSpectra"))
    res_2 <- lapply(spl, removePeaks)
    expect_equal(res@listData, res_2)
    expect_equal(mcols(res)$id, 1:3)

    res <- removePeaks(spl, t = 10, msLevel. = 1)
    expect_true(is(res, "MSpectra"))
    res_2 <- lapply(spl, removePeaks, t = 10, msLevel. = 1)
    expect_equal(res@listData, res_2)
    expect_equal(res[[3]], removePeaks(spl[[3]], t = 10))

    msnexp_remp <- removePeaks(tmt_od_sub)
    res_spctra <- removePeaks(MSpectra(spectra(tmt_od_sub)))
    expect_equal(spectra(msnexp_remp), res_spctra@listData)
})

test_that("filterMz,MSpectra works", {
    spl <- MSpectra(new("Spectrum2", mz = 1:5, intensity = c(4, 34, 5, 199, 3),
                       centroided = TRUE),
                   new("Spectrum2", mz = 1:4, intensity = c(123, 43, 3, 9),
                       centroided = TRUE),
                   new("Spectrum1", mz = 1:4, intensity = c(123, 43, 2, 9),
                       centroided = TRUE),
                   elementMetadata = DataFrame(id = 1:3))
    res <- filterMz(spl, mz = c(3, 5))
    expect_true(all(unlist(mz(res)) >= 3))
    expect_equal(res@listData, lapply(spl, filterMz, mz = c(3, 5)))
    res <- filterMz(spl, mz = c(3, 5), msLevel. = 1)
    expect_equal(res[1], spl[1])
    expect_equal(res[[3]], filterMz(spl[[3]], mz = c(3, 5)))

    msnexp_filt <- filterMz(tmt_od_sub, mz = c(500, 600))
    expect_warning(res_spctra <- filterMz(MSpectra(spectra(tmt_od_sub)),
                                          mz = c(500, 600)))
    expect_equal(spectra(msnexp_filt), res_spctra@listData)
})

test_that("pickPeaks,MSpectra and smooth,MSpectra works", {
    spctra <- spectra(filterFile(sciex, 1))
    spl <- MSpectra(spctra)
    res <- pickPeaks(spl)
    expect_true(is(res, "MSpectra"))
    expect_equal(res@listData, lapply(spctra, pickPeaks))
    expect_true(all(peaksCount(res) < peaksCount(spl)))

    res <- pickPeaks(spl, msLevel = 2:3)
    expect_equal(peaksCount(res), peaksCount(spl))

    res <- pickPeaks(spl, msLevel = 1:4)
    expect_true(all(peaksCount(res) < peaksCount(spl)))

    expect_warning(res <- smooth(spl))
    expect_true(is(res, "MSpectra"))
    expect_equal(res@listData, lapply(spctra, smooth))
})

test_that("combineSpectra,MSpectra works", {
    set.seed(123)
    mzs <- seq(1, 20, 0.1)
    ints1 <- abs(rnorm(length(mzs), 10))
    ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
    ints2 <- abs(rnorm(length(mzs), 10))
    ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
    ints3 <- abs(rnorm(length(mzs), 10))
    ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)

    ## Create the spectra.
    sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
               intensity = ints1, rt = 1)
    sp2 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
               intensity = ints2, rt = 2)
    sp3 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.009),
               intensity = ints3, rt = 3)
    sp4 <- new("Spectrum2", mz = 1:3, intensity = c(4, 8, 1))
    sp5 <- new("Spectrum2", mz = 2:4, intensity = c(5, 8, 2))

    spctra <- MSpectra(sp1, sp2, sp3,
                      elementMetadata = DataFrame(idx = 1:3,
                                                  group = c("b", "a", "a")))

    res <- combineSpectra(spctra)
    expect_true(length(res) == 1)
    expect_equal(rtime(res), c(`1` = 1))
    expect_equal(mcols(res), DataFrame(idx = 1, group = "b", row.names = "1"))

    names(spctra) <- c("A", "B", "C")
    res <- combineSpectra(spctra)

    expect_error(combineSpectra(spctra, fcol = "other"))
    res <- combineSpectra(spctra, fcol = "group", mzd = 0.05)
    expect_equal(lengths(intensity(res)), c(A = 191, B = 191))
    expect_equal(names(res), c("A", "B"))
    expect_equal(res[[1]], sp1)
    expect_equal(mcols(res), DataFrame(idx = c(1, 2), group = c("b", "a"),
                                       row.names = c("A", "B")))

    spctra <- MSpectra(sp1, sp2, sp3, sp4, sp5,
                      elementMetadata = DataFrame(group = c("a", "b", "a",
                                                            "c", "c")))
    expect_error(combineSpectra(spctra))
    res <- combineSpectra(spctra, fcol = "group", mzd = 0.05)
    expect_equal(names(res), c("1", "2", "4"))
    expect_equal(res[[2]], sp2)
    expect_equal(rtime(res), c(`1` = 1, `2` = 2, `4` = NA))
    expect_equal(msLevel(res), c(`1` = 1, `2` = 1, `4` = 2))
    expect_equal(intensity(res)[[3]], c(4, mean(c(8, 5)), mean(c(1, 8)), 2))
})

test_that("as,MSpectra,list works", {
    sps <- MSpectra()
    res <- as(sps, "list")
    expect_true(length(res) == 0)
    expect_true(is.list(res))

    res <- as(spl_, "list")
    expect_true(is.list(res))
    expect_equivalent(res, list(sp1, sp2, sp3))
})

test_that("as,MSpectra,MSnExp works", {
    sps <- MSpectra()
    res <- as(sps, "MSnExp")
    expect_true(is(res, "MSnExp"))
    expect_true(length(res) == 0)
    expect_true(validObject(res))

    expect_error(as(spl_, "MSnExp"))
    spl_2 <- filterMsLevel(spl_, 2)
    res <- as(spl_2, "MSnExp")
    expect_equal(spectra(res), as(spl_2, "list"))

    ## res <- as(spl_2[2:1], "MSnExp")
    ## expect_equal(spectra(res), as(spl_2[2:1], "list"))
    ## mcols(spl_2) <- DataFrame(id = c("a", "b"), other = c("d", "e"))

    res <- as(spl_2, "MSnExp")
    expect_equal(fData(res)$id, mcols(spl_2)$id)
    expect_equal(fData(res)$other, mcols(spl_2)$other)
    expect_equal(spectra(res), as(spl_2, "list"))

    spl_1 <- filterMsLevel(spl_, 1)
    res <- as(spl_1, "MSnExp")
    expect_equal(spectra(res), list(`3` = sp3))
})

test_that("filterMsLevel,MSpectra works", {
    res <- filterMsLevel(spl_, 1)
    sp_3 <- MSpectra(sp3)
    names(sp_3) <- "3"
    expect_equal(res, sp_3)

    expect_equal(filterMsLevel(spl_), spl_)
    expect_equal(filterMsLevel(spl_, 1:4), spl_)

    res <- filterMsLevel(spl_, 2)
    expect_equal(res, MSpectra(sp1, sp2))

    res <- filterMsLevel(spl_, 4)
    expect_true(is(res, "MSpectra"))
    expect_true(length(res) == 0)
})
