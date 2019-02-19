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

test_that(".valid.MSnExperiment.featureData works", {
    vld <- DataFrame(fileIdx = c(1L, 1L, 2L, 2L, 3L),
                     msLevel = c(1L, 1L, 2L, 2L, 1L))
    expect_null(.valid.MSnExperiment.featureData(vld))
    expect_match(.valid.MSnExperiment.featureData(vld, 2), "number of distinct")
    vld$msLevel[4] <- NA
    expect_match(.valid.MSnExperiment.featureData(vld),
                 "'msLevel' should contain")
    vld$msLevel[4] <- 4.23
    expect_match(.valid.MSnExperiment.featureData(vld),
                 "'msLevel' should contain")
    vld$msLevel[4] <- 4L
    vld$msLevel <- as.integer(vld$msLevel)
    vld$fileIdx[1] <- NA
    expect_match(.valid.MSnExperiment.featureData(vld), "'fileIdx' should")
})

test_that("MSnExperiment constructor works", {
    ## From a list of spectra.
    sp1 <- new("Spectrum1", rt = 1.2, mz = 1:4, intensity = abs(rnorm(4)))
    sp2 <- new("Spectrum1", rt = 1.4, mz = 1:4, intensity = abs(rnorm(4)))
    sp3 <- new("Spectrum1", rt = 1.5, mz = 1:4, intensity = abs(rnorm(4)))
    spl <- list(sp1, sp2, sp3)
    res <- MSnExperiment(spl)
    expect_equal(rownames(spectraData(res)), c("F1.S1", "F1.S2", "F1.S3"))
    sp1@fromFile <- 1L
    sp2@fromFile <- 1L
    sp3@fromFile <- 1L
    expect_equal(res[[1]], sp1)
    expect_equal(res[[2]], sp2)
    expect_equal(res[[3]], sp3)

    ## Providing a sampleData
    res <- MSnExperiment(spl, sampleData = DataFrame(a = "some file"))
    expect_equal(res@sampleData, DataFrame(a = "some file"))

    ## Providing spectraData
    spd <- DataFrame(peak_id = c("CP1", "CP1", "CP2"),
                     compound_id = c("HMBD1", "HMDB1", "HMDB2"))
    res <- MSnExperiment(spl, spectraData = spd)
    expect_equal(res@spectraData$peak_id, spd$peak_id)
    expect_equal(res@spectraData$compound_id, spd$compound_id)

    ## Using named list
    names(spl) <- c("A", "B", "C")
    res <- MSnExperiment(spl)
    expect_equal(rownames(res@spectraData), c("A", "B", "C"))

    ## Errors
    expect_error(MSnExperiment(4), "'x' has to be")
    expect_error(MSnExperiment(spl, spectraData = DataFrame(iam = "wrong")))
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

test_that("removePeaks,MSnExperiment and clean,MSnExperiment work", {
    sciex_spctra <- spectra(sciex)
    tmp <- removePeaks(sciex_mzr, t = 10000)
    tmp_inmem <- removePeaks(sciex_inmem, t = 10000)
    expect_true(length(tmp@processingQueue) == 1)
    sciex_spctra <- lapply(sciex_spctra, removePeaks, t = 10000)
    expect_equal(spectrapply(tmp), sciex_spctra)
    expect_equal(spectrapply(tmp_inmem), sciex_spctra)

    tmp <- clean(tmp, all = TRUE)
    tmp_inmem <- clean(tmp_inmem, all = TRUE)
    expect_true(length(tmp@processingQueue) == 2)
    sciex_spctra <- lapply(sciex_spctra, clean, all = TRUE)
    expect_equal(spectrapply(tmp), sciex_spctra)
    expect_equal(spectrapply(tmp_inmem), sciex_spctra)
})

test_that("fileNames,MSnExperiment works", {
    expect_equal(unname(fileNames(sciex_mzr)), sf)
    expect_equal(fileNames(new("MSnExperiment", backend = BackendMzR())),
                 character())
})

test_that("filterFile and [,MSnExperiment work", {
    fls <-  c(system.file("microtofq/MM14.mzML", package = "msdata"),
              system.file("microtofq/MM8.mzML", package = "msdata"),
              sf)
    mse <- readMSnExperiment(fls, backend = BackendMzR())
    mse_hdf5 <- readMSnExperiment(fls, backend = BackendHdf5(),
                                  path = paste0(tempdir(), "/mse/"))
    mse_mem <- readMSnExperiment(fls, backend = BackendMemory())
    sps <- spectrapply(mse_mem)
    ## filterFile
    res <- filterFile(mse, c(FALSE, TRUE, TRUE, FALSE))
    res_hdf5 <- filterFile(mse_hdf5, c(FALSE, TRUE, TRUE, FALSE))
    res_mem <- filterFile(mse_mem, c(FALSE, TRUE, TRUE, FALSE))
    expect_equal(unname(fileNames(res)), fls[2:3])
    expect_equal(fileNames(res), fileNames(res_hdf5))
    expect_equal(fileNames(res), fileNames(res_mem))
    expect_equal(res@sampleData, mse@sampleData[2:3, , drop = FALSE])
    expect_equal(res@sampleData, res_hdf5@sampleData)
    expect_equal(res@sampleData, res_mem@sampleData)
    expect_equal(res[[1]]@fromFile, 1L)
    expect_equal(res_hdf5[[1]]@fromFile, 1L)
    expect_equal(res_mem[[1]]@fromFile, 1L)
    expect_equal(spectrapply(res, FUN = intensity),
                 lapply(sps[mse@spectraData$fileIdx %in% 2:3], intensity))
    expect_equal(spectrapply(res), spectrapply(res_hdf5))
    expect_equal(spectrapply(res), spectrapply(res_mem))

    ## filterFile, change order
    res <- filterFile(mse, c(2, 1))
    res_hdf5 <- filterFile(mse_hdf5, c(2, 1))
    res_mem <- filterFile(mse_mem, c(2, 1))
    expect_equal(unname(fileNames(res)), fls[2:1])
    expect_equal(fileNames(res), fileNames(res_hdf5))
    expect_equal(fileNames(res), fileNames(res_mem))
    expect_equal(res@sampleData, mse@sampleData[2:1, , drop = FALSE])
    expect_equal(res@sampleData, res_hdf5@sampleData)
    expect_equal(res@sampleData, res_mem@sampleData)
    expect_equal(spectrapply(res, FUN = intensity),
                 lapply(sps[c(113:310, 1:112)], intensity))
    expect_equal(spectrapply(res), spectrapply(res_hdf5))
    expect_equal(spectrapply(res), spectrapply(res_mem))

    ## [, subset spectra - in arbitrary order.
    idx <- c(333, 323, 17, 21, 337) # file 3, file 3, file 1, file 1, file 3
    res <- mse[idx]
    res_hdf5 <- mse_hdf5[idx]
    res_mem <- mse_mem[idx]
    expect_equal(unname(fileNames(res)), fls[c(3, 1)])
    expect_equal(fileNames(res), fileNames(res_hdf5))
    expect_equal(fileNames(res), fileNames(res_mem))
    expect_equal(res@spectraData$fileIdx, c(1, 1, 2, 2, 1))
    expect_equal(res@sampleData, mse@sampleData[c(3, 1), , drop = FALSE])
    expect_equal(res@sampleData, res_mem@sampleData)
    expect_equal(res@sampleData, res_hdf5@sampleData)
    res_sps <- spectrapply(res)
    expect_equal(names(res_sps), names(sps)[idx])
    expect_equal(lapply(res_sps, intensity),
                 lapply(sps[idx], intensity))
    expect_equal(res_sps, spectrapply(res_hdf5))
    expect_equal(res_sps, spectrapply(res_mem))
    ## drop
    res <- mse[333]
    expect_true(is(res, "Spectrum"))
    expect_equal(fromFile(res), 1)
    expect_equal(intensity(res), intensity(sps[[333]]))
    expect_equal(res, mse_hdf5[333])
    expect_equal(res, mse_mem[333])
    res <- mse[333, drop = FALSE]
    expect_true(is(res, "MSnExperiment"))

    ## Errors
    expect_error(mse[1, 4])
    expect_error(mse[, 5])
    expect_error(mse[12222222, ])
})

test_that("[[,MSnExperiment works", {
    res <- sciex_mzr[[13]]
    expect_true(inherits(res, "Spectrum"))
    expect_equal(sciex[[13]], res)
    expect_error(sciex_mzr[[1222222]])
})

test_that("filterFile works", {
    res <- filterFile(sciex_mzr, 2)
    expect_equal(fileNames(res), fileNames(sciex_mzr)[2])
    expect_true(all(res@spectraData$fileIdx == 1))
    expect_equal(res@sampleData, sciex_mzr@sampleData[2, , drop = FALSE])
    expect_equal(spectrapply(res, return.type = "list"),
                 spectrapply(filterFile(sciex, 2)))
    expect_error(filterFile(sciex_mzr, "b"))
    expect_error(filterFile(sciex_mzr, c(1, 1, 1, 2)))
    expect_error(filterFile(sciex_mzr, 5))
})

test_that("setBackend methods work", {
    ## MzR -> Memory
    sciex_mzr_mem <- setBackend(sciex_mzr, backend = BackendMemory())
    expect_true(is(sciex_mzr_mem@backend, "BackendMemory"))
    expect_true(validObject(sciex_mzr_mem@backend))
    expect_equal(spectrapply(sciex_mzr_mem), spectrapply(sciex_inmem))

    ## MzR -> Hdf5
    sciex_mzr_h5 <- setBackend(sciex_mzr, backend = BackendHdf5(),
                                  path = paste0(tempdir(), "/switch1/"))
    expect_true(is(sciex_mzr_h5@backend, "BackendHdf5"))
    expect_true(validObject(sciex_mzr_h5@backend))
    expect_equal(spectrapply(sciex_mzr_h5), spectrapply(sciex_inmem))
    ## Memory -> MzR
    sciex_mem_mzr <- setBackend(sciex_inmem, backend = BackendMzR())
    expect_true(is(sciex_mem_mzr@backend, "BackendMzR"))
    expect_true(validObject(sciex_mem_mzr@backend))
    expect_equal(spectrapply(sciex_mem_mzr), spectrapply(sciex_inmem))

    ## Memory, modify -> MzR
    tmp <- removePeaks(sciex_inmem, t = 5000)

    ## Memory, modify -> Hdf5

    ## The tests below are redundant
    ## ## Memory -> Hdf5
    ## sciex_mem_h5 <- setBackend(sciex_inmem, backend = BackendHdf5(),
    ##                               path = paste0(tempdir(), "/switch2/"))
    ## expect_true(is(sciex_mem_h5@backend, "BackendHdf5"))
    ## expect_true(validObject(sciex_mem_h5@backend))
    ## expect_equal(spectra(sciex_mem_h5), spectra(sciex_inmem))
    ## ## Hdf5 -> Memory
    ## sciex_h5_mem <- setBackend(sciex_h5, backend = BackendMemory())
    ## expect_true(is(sciex_h5_mem@backend, "BackendMemory"))
    ## expect_true(validObject(sciex_h5_mem@backend))
    ## expect_equal(spectra(sciex_h5_mem), spectra(sciex_inmem))
    ## ## Hdf5 -> MzR
    ## sciex_h5_mzr <- setBackend(sciex_h5, backend = BackendMzR())
    ## expect_true(is(sciex_h5_mzr@backend, "BackendMzR"))
    ## expect_true(validObject(sciex_h5_mzr@backend))
    ## expect_equal(spectra(sciex_h5_mzr), spectra(sciex_inmem))
})

test_that("as,MSnExperiment works", {
    subs <- sciex_inmem[1:34]
    res <- as(subs, "list")
    expect_true(is.list(res))
    expect_equal(res, spectrapply(subs))
    res <- as(subs, "List")
    expect_true(is(res, "List"))

    res <- as(sciex_h5[1:34], "list")
    expect_equal(res, as(subs, "list"))

    res <- as(sciex_mzr[1:34], "list")
    expect_equal(res, as(subs, "list"))
})

test_that("featureData, and featureData<-,MSnExperiment work", {
    spl <- list(new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))),
                new("Spectrum1", mz = 1:4, intensity = abs(rnorm(4))),
                new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))))
    mse <- MSnExperiment(spl)
    featureData(mse)$polarity <- c(-1, 1, 1)
    expect_equal(featureData(mse)$polarity, c(-1, 1, 1))
    expect_equal(unname(vapply(as(mse, "list"), polarity, integer(1))),
                 c(-1, 1, 1))
    expect_error(featureData(mse)$fileIdx <- c("a", "b", "c"))
    expect_error(featureData(mse) <- DataFrame(a = 1:3))

    featureData(mse)$feature_id <- c("FT01", "FT01", "FT02")
    expect_equal(featureData(mse)$feature_id, c("FT01", "FT01", "FT02"))
    expect_equal(unname(spectrapply(mse, FUN = intensity)),
                 lapply(spl, intensity))

    spectraData(mse)$peak_id <- 1:3
    expect_equal(spectraData(mse)$peak_id, 1:3)
})

test_that("sampleData, sampleData<-,MSnExperiment work", {
    spl <- list(new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))),
                new("Spectrum1", mz = 1:4, intensity = abs(rnorm(4))),
                new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))))
    mse <- MSnExperiment(spl)
    sampleData(mse) <- DataFrame(sample_name = "a", sample_id = "001")
    expect_equal(sampleData(mse)$sample_name, "a")
    expect_equal(sampleData(mse)$sample_id, "001")
    expect_error(sampleData(mse) <- "a")
    expect_error(sampleData(mse) <- DataFrame(sample_name = c("a", "b")))
})

test_that("length,MSnExperiment works", {
    expect_equal(length(sciex_mzr), 1862)
})

test_that("featureNames,MSnExperiment works", {
    expect_equal(featureNames(sciex_mzr), rownames(sciex_mzr@spectraData))
})

test_that("acquisitionNum,MSnExperiment works", {
    spl <- list(new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))),
                new("Spectrum1", mz = 1:4, intensity = abs(rnorm(4))),
                new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))))
    mse <- MSnExperiment(spl)
    expect_equal(acquisitionNum(mse), c(F1.S1 = NA_integer_,
                                        F1.S2 = NA_integer_,
                                        F1.S3 = NA_integer_))
    expect_equal(unname(acquisitionNum(sciex_mzr)),
                 sciex_mzr@spectraData$acquisitionNum)
    expect_equal(names(acquisitionNum(sciex_mzr)),
                 rownames(sciex_mzr@spectraData))
})
