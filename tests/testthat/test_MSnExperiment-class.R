test_that("MSnExperiment validator works", {
    tst <- new("MSnExperiment")
    expect_error(validObject(tst))
    tst <- new("MSnExperiment", backend = BackendMzR())
    expect_true(validObject(tst))
    tst@sampleData <- DataFrame(sampleIdx = 1)
    res <- MSnbase:::.validMSnExperiment(tst)
    expect_true(is.character(res))
    expect_error(validObject(tst))
})

test_that("addProcessingStep works", {
    tst <- new("MSnExperiment", backend = BackendMzR())
    tst <- addProcessingStep(tst, mean)
    expect_true(length(tst@processingQueue) == 1)
    expect_error(addProcessingStep(tst, "4"))
    tst <- addProcessingStep(tst, function(z, t) z * t, t = 4)
    expect_true(length(tst@processingQueue) == 2)
})

test_that(".valid.processingQueue works", {
    expect_true(is.character(.valid.processingQueue(list(3, 5))))
    lst <- list(ProcessingStep(mean), ProcessingStep("max"))
    expect_true(is.null(.valid.processingQueue(lst)))
})

test_that("readMSnExperiment works", {
    sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    res <- readMSnExperiment(sf)
    expect_true(validObject(res))
    expect_true(length(unique(res@spectraData$fileIdx)) == 2)
    expect_error(readMSnExperiment(sf, sampleData = DataFrame(a = 3)))
})

test_that("spectraData works", {
    expect_error(spectraData(4))
    expect_equal(spectraData(sciex_mzr), sciex_mzr@spectraData)
})

test_that("spectrapply,MSnExperiment works", {
    sciex_spctra <- spectra(sciex)

    expect_error(spectrapply(sciex_mzr, f = 4))

    ## BackendMzR
    res <- spectrapply(sciex_mzr)
    expect_equal(names(res), rownames(sciex_mzr@spectraData))
    expect_true(all(vapply(res, is, logical(1), "Spectrum")))
    expect_true(all(vapply(res, validObject, logical(1))))
    expect_equal(sciex_spctra, res)

    ## No parallelization:
    expect_equal(sciex_spctra,
                 spectrapply(sciex_mzr, f = rep(1, nrow(sciex_mzr@spectraData))))

    res_2 <- spectrapply(sciex_mzr, FUN = function(z) mean(mz(z)))
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))

    sciex_mzr@processingQueue <- list(ProcessingStep("mz"))
    res_2 <- spectrapply(sciex_mzr, FUN = "mean")
    expect_true(all(vapply(res_2, is.numeric, logical(1))))
    expect_equal(res_2, lapply(sciex_spctra, function(z) mean(mz(z))))
    sciex_mzr@processingQueue <- list()

    ## Work on subsets - to simulate subsetting/filtering
    idx <- c(4, 123, 432, 1110, 1234)
    tmp <- sciex_mzr
    tmp@spectraData <- tmp@spectraData[idx, ]
    res <- spectrapply(tmp)
    expect_equal(res, sciex_spctra[idx])

    ## Compare BackendMzR and BackendMemory
    expect_equal(spectrapply(sciex_inmem), spectrapply(sciex_mzr))
    expect_equal(spectrapply(sciex_inmem, FUN = mz),
                 spectrapply(sciex_mzr, FUN = mz))
})

test_that("spectra,MSnExperiment works", {
    sciex_spctra <- spectra(sciex)
    expect_equal(spectra(sciex_mzr, return.type = "list"), sciex_spctra)
})

test_that("removePeaks,MSnExperiment and clean,MSnExperiment work", {
    sciex_spctra <- spectra(sciex)
    tmp <- removePeaks(sciex_mzr, t = 10000)
    tmp_inmem <- removePeaks(sciex_inmem, t = 10000)
    expect_true(length(tmp@processingQueue) == 1)
    sciex_spctra <- lapply(sciex_spctra, removePeaks, t = 10000)
    expect_equal(spectra(tmp, return.type = "list"), sciex_spctra)
    expect_equal(spectra(tmp_inmem, return.type = "list"), sciex_spctra)

    tmp <- clean(tmp, all = TRUE)
    tmp_inmem <- clean(tmp_inmem, all = TRUE)
    expect_true(length(tmp@processingQueue) == 2)
    sciex_spctra <- lapply(sciex_spctra, clean, all = TRUE)
    expect_equal(spectra(tmp, return.type = "list"), sciex_spctra)
    expect_equal(spectra(tmp_inmem, return.type = "list"), sciex_spctra)
})

test_that("fileNames,MSnExperiment works", {
    expect_equal(fileNames(sciex_mzr), sf)
    expect_equal(fileNames(new("MSnExperiment", backend = BackendMzR())),
                 character())
})

test_that("[,MSnExperiment works", {
    fls <-  c(system.file("microtofq/MM14.mzML", package = "msdata"),
              system.file("microtofq/MM8.mzML", package = "msdata"),
              sf)
    mse <- readMSnExperiment(fls)
    spctra <- spectra(mse)
    ## Subset files.
    res <- mse[, c(FALSE, TRUE, TRUE, FALSE)]
    expect_equal(fileNames(res), fls[2:3])
    expect_equal(spectrapply(res, FUN = mz),
                 lapply(spctra[mse@spectraData$fileIdx %in% 2:3], mz))
    expect_equal(res@sampleData, mse@sampleData[2:3, , drop = FALSE])
    ## Subset files, change order
    res <- mse[, c(2, 1)]
    expect_equal(fileNames(res), fls[2:1])
    expect_equal(res@sampleData, mse@sampleData[c(2, 1), , drop = FALSE])
    expect_equal(spectrapply(res, FUN = intensity),
                 lapply(spctra[c(113:310, 1:112)], intensity))
    ## Subset spectra - in arbitrary order.
    idx <- c(333, 323, 17, 21, 337) # file 3, file 3, file 2, file 2, file 3
    res <- mse[idx, ]
    expect_equal(fileNames(res), fls[c(3, 1)])
    expect_equal(res@spectraData$fileIdx, c(1, 1, 2, 2, 1))
    res_spctra <- spectra(res)
    expect_equal(names(res_spctra), names(spctra)[idx])
    expect_equal(lapply(res_spctra, intensity),
                 lapply(spctra[idx], intensity))
    ## drop
    res <- mse[333]
    expect_true(is(res, "Spectrum"))
    expect_equal(fromFile(res), 1)
    expect_equal(intensity(res), intensity(spctra[[333]]))
    res <- mse[333, , drop = FALSE]
    expect_true(is(res, "MSnExperiment"))

    ## Errors
    expect_error(mse[1, 4])
    expect_error(mse[, 5])
    expect_error(mse[12222222, ])
})
