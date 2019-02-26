test_that("MSnExperiment validator works", {
    tst <- new("MSnExperiment")
    expect_error(validObject(tst))
    tst <- new("MSnExperiment", backend = BackendMzR())
    expect_true(validObject(tst))
    tst@sampleData <- DataFrame(sampleIdx = 1)
    res <- .validMSnExperiment(tst)
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

    ## From empty spectra.
    mse <- MSnExperiment(list(new("Spectrum1"), new("Spectrum1", mz = 1:3,
                                                    intensity = 1:3)))
    mse <- MSnExperiment(list(new("Spectrum2"), new("Spectrum1", mz = 1:3,
                                                    intensity = 1:3)))

    ## From multiple files
    mse <- MSnExperiment(list(new("Spectrum1", fromFile = 2L),
                              new("Spectrum1", fromFile = 1L),
                              new("Spectrum1", fromFile = 1L),
                              new("Spectrum1", fromFile = 2L),
                              new("Spectrum1", fromFile = 2L)))
    expect_equal(nrow(sampleData(mse)), 2)
    expect_equal(featureNames(mse), c("F2.S1", "F1.S1", "F1.S2", "F2.S2",
                                      "F2.S3"))

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
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
           system.file("microtofq/MM8.mzML", package = "msdata"))
    mzr <- readMSnExperiment(f, backend = BackendMzR())
    mem <- readMSnExperiment(f, backend = BackendMemory())

    ## MzR -> Memory
    mzr_mem <- setBackend(mzr, backend = BackendMemory())
    expect_true(is(mzr_mem@backend, "BackendMemory"))
    expect_true(validObject(mzr_mem@backend))
    expect_equal(spectrapply(mzr_mem), spectrapply(mem))
    ## MzR -> Hdf5
    mzr_h5 <- setBackend(mzr, backend = BackendHdf5(),
                         path = paste0(tempdir(), "/switch1/"))
    expect_true(is(mzr_h5@backend, "BackendHdf5"))
    expect_true(validObject(mzr_h5@backend))
    expect_equal(spectrapply(mzr_h5), spectrapply(mem))
    ## Memory -> MzR
    mem_mzr <- setBackend(mem, backend = BackendMzR())
    expect_true(is(mem_mzr@backend, "BackendMzR"))
    expect_true(validObject(mem_mzr@backend))
    expect_equal(spectrapply(mem_mzr), spectrapply(mem))

    ## Memory, modify -> MzR
    tmp <- removePeaks(mem, t = 5000)
    tmp <- applyProcessingQueue(tmp)
    expect_equal(tmp@processingQueue, list())
    expect_equal(tmp@backend@modCount, c(1L, 1L))
    expect_error(setBackend(tmp, BackendMzR()), "Can not change backend")

    ## Memory, modify -> Hdf5
    tmp_h5 <- setBackend(tmp, BackendHdf5(),
                         path = paste0(tempdir(), "/switch2/"))
    expect_true(length(tmp_h5@processingQueue) == 0)
    expect_equal(unname(tmp_h5@backend@modCount), c(1L, 1L))
    sps_h5 <- spectrapply(tmp_h5)
    sps <- spectrapply(mem)
    expect_equal(lapply(sps_h5, ionCount),
                 lapply(sps, function(z) ionCount(removePeaks(z, t = 5000))))
    expect_error(setBackend(tmp_h5, BackendMzR()), "Can not change backend")

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

test_that("applyProcessingQueue works", {
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
           system.file("microtofq/MM8.mzML", package = "msdata"))
    ## Memory backend
    mse_mem <- readMSnExperiment(f, backend = BackendMemory())
    sps <- spectrapply(mse_mem)
    idx <- c(13, 15, 33, 89, 113, 117, 167)
    mse_mem <- mse_mem[idx]
    mse_mem <- removePeaks(mse_mem, t = 5000)
    expect_equal(mse_mem@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_mem@processingQueue), 1L)
    mse_mem <- applyProcessingQueue(mse_mem)
    expect_equal(mse_mem@backend@modCount, c(1L, 1L))
    expect_true(length(mse_mem@processingQueue) == 0)
    expect_equal(as(mse_mem, "list"), lapply(sps[idx], removePeaks, t = 5000))
    ## Hdf5 backend
    mse_h5 <- readMSnExperiment(f, backend = BackendHdf5(),
                                path = paste0(tempdir(), "/apq"))
    mse_h5 <- mse_h5[idx]
    sps <- spectrapply(mse_h5)
    expect_equal(mse_h5@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_h5@processingQueue), 0L)
    mse_h5 <- removePeaks(mse_h5, t = 5000)
    expect_equal(mse_h5@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_h5@processingQueue), 1L)
    mse_h5 <- applyProcessingQueue(mse_h5)
    expect_equal(mse_h5@backend@modCount, c(1L, 1L))
    expect_equal(length(mse_h5@processingQueue), 0L)
    expect_equal(spectrapply(mse_h5, ionCount),
                 lapply(sps, function(z) ionCount(removePeaks(z, t = 5000))))
    ## MzR backend
    mse_mzr <- readMSnExperiment(f, backend = BackendMzR())
    mse_mzr <- mse_mzr[idx]
    sps <- spectrapply(mse_mzr)
    expect_equal(mse_mzr@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_mzr@processingQueue), 0L)
    mse_mzr <- removePeaks(mse_mzr, t = 5000)
    expect_equal(mse_mzr@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_mzr@processingQueue), 1L)
    mse_mzr <- applyProcessingQueue(mse_mzr)
    expect_equal(mse_mzr@backend@modCount, c(0L, 0L))
    expect_equal(length(mse_mzr@processingQueue), 1L)
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

test_that("featureNames, spectraNames,MSnExperiment work", {
    expect_equal(featureNames(sciex_mzr), rownames(sciex_mzr@spectraData))
    expect_equal(spectraNames(sciex_mzr), rownames(sciex_mzr@spectraData))
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

test_that("centroided, isCentroided work", {
    spl <- list(new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))),
                new("Spectrum1", mz = 1:4, intensity = abs(rnorm(4))),
                new("Spectrum1", mz = 1:5, intensity = abs(rnorm(5))))
    mse <- MSnExperiment(spl)
    expect_true(all(is.na(centroided(mse))))
    expect_equal(names(centroided(mse)), rownames(mse@spectraData))
    expect_error(centroided(mse, na.fail = TRUE))
    expect_error(centroided(mse) <- "a")
    expect_error(centroided(mse) <- c(TRUE, FALSE))
    centroided(mse) <- FALSE
    expect_true(all(!centroided(mse)))
    centroided(mse) <- c(TRUE, FALSE, TRUE)
    expect_equal(unname(centroided(mse)), c(TRUE, FALSE, TRUE))

    centroided(sciex_mzr) <- TRUE
    expect_true(centroided(sciex_mzr[[13]]))
    centroided(sciex_mzr) <- FALSE
    tmp <- sciex_inmem[1:10]
    centroided(tmp) <- TRUE
    expect_true(centroided(tmp[[2]]))

    expect_true(all(!isCentroided(tmp)))
})

test_that("collisionEnergy,collisionEnergy<-,MSnExperiment work", {
    tmp <- sciex_inmem[1:10]
    expect_true(all(collisionEnergy(tmp) == 0))
    expect_equal(names(collisionEnergy(tmp)), rownames(tmp@spectraData))

    expect_error(collisionEnergy(tmp) <- c(1.2, 1.5))
    expect_error(collisionEnergy(tmp) <- rep("a", 10))
    collisionEnergy(tmp) <- 1:10
})

test_that("fromFile,MSnExperiment works", {
    tmp <- sciex_inmem[1:10]
    expect_true(all(fromFile(tmp) == 1))
    expect_equal(names(fromFile(tmp)), rownames(tmp@spectraData))
})

test_that("intensity and mz,MSnExperiment work", {
    tmp <- sciex_inmem[1:4]
    res <- intensity(tmp)
    expect_equal(names(res), rownames(tmp@spectraData))
    expect_true(all(vapply(res, is.numeric, logical(1))))
    expect_equal(unname(lengths(res)), c(578, 1529, 1600, 1664))
    res <- mz(tmp)
    expect_equal(names(res), rownames(tmp@spectraData))
    expect_true(all(vapply(res, is.numeric, logical(1))))
    expect_false(all(vapply(res, is.unsorted, logical(1))))
    expect_equal(unname(lengths(res)), c(578, 1529, 1600, 1664))

    mse <- MSnExperiment(list(new("Spectrum2"), new("Spectrum2")))
    res <- intensity(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_true(all(lengths(res) == 0))
    res <- mz(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_true(all(lengths(res) == 0))
})

test_that("ionCount, and isEmpty,MSnExperiment work", {
    tmp <- sciex_inmem[1:4]
    res <- ionCount(tmp)
    expect_equal(names(res), rownames(tmp@spectraData))
    expect_true(is.numeric(res))
    tmp <- removePeaks(tmp, t = 5000)
    res_2 <- ionCount(tmp)
    expect_true(all(res_2 < res))

    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2"),
                              new("Spectrum1", mz = 1:3, intensity = 1:3)))
    res <- ionCount(mse)
    expect_equal(unname(res), c(10, 0, 6))

    res <- isEmpty(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_true(is.logical(res))
    expect_equal(unname(res), c(FALSE, TRUE, FALSE))
})

test_that("msLevel,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2"),
                              new("Spectrum1", mz = 1:3, intensity = 1:3)))
    expect_equal(unname(msLevel(mse)), c(1L, 2L, 1L))
    expect_equal(names(msLevel(mse)), rownames(mse@spectraData))
})

test_that("polarity, polarity<-,MSnExperiment work", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2"),
                              new("Spectrum1", mz = 1:3, intensity = 1:3)))
    expect_equal(unname(polarity(mse)), rep(NA_integer_, 3))
    expect_equal(names(msLevel(mse)), rownames(mse@spectraData))
    polarity(mse) <- 1L
    expect_equal(unname(polarity(mse)), rep(1L, 3))
    polarity(mse) <- c(1L, 0L, 1L)
    expect_equal(unname(polarity(mse)), c(1L, 0L, 1L))
    expect_equal(unname(polarity(spectrapply(mse)[[3]])), 1L)
    expect_error(polarity(mse) <- "a")
    expect_error(polarity(mse) <- c(2L, 1L))
})

test_that("rtime, rtime<- MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2"),
                              new("Spectrum1", mz = 1:3, intensity = 1:3)))
    expect_equal(unname(rtime(mse)), rep(NA_real_, 3))
    expect_equal(names(rtime(mse)), rownames(mse@spectraData))
    rtime(mse) <- c(1.2, 1.3, 1.4)
    expect_equal(unname(rtime(mse)), c(1.2, 1.3, 1.4))
    expect_equal(unname(rtime(spectrapply(mse)[[1]])), 1.2)
    expect_error(rtime(mse) <- "f")
    expect_error(rtime(mse) <- c(3, 4))
})

test_that("peaksCount,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2"),
                              new("Spectrum1", mz = 1:3, intensity = 1:3)))
    expect_equal(unname(peaksCount(mse)), c(4L, 0, 3L))
    expect_equal(names(peaksCount(mse)), rownames(mse@spectraData))
})

test_that("precursor*,MSnExperiment methods work", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4),
                              new("Spectrum2", precScanNum = 1L,
                                  precursorIntensity = 3, precursorMz = 3),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  precursorCharge = 1L, precursorIntensity = 2,
                                  precScanNum = 1L)))
    res <- precursorCharge(mse)
    expect_equal(res, c(F1.S1 = 0L, F1.S2 = NA_integer_, F1.S3 = 1L))
    res <- precursorIntensity(mse)
    expect_equal(res, c(F1.S1 = 0, F1.S2 = 3, F1.S3 = 2))
    res <- precursorMz(mse)
    expect_equal(res, c(F1.S1 = 0, F1.S2 = 3, F1.S3 = NA))
    res <- precScanNum(mse)
    expect_equal(res, c(F1.S1 = 0L, F1.S2 = 1L, F1.S3 = 1L))
})

test_that("bpi, tic, MSnExperiment work", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  tic = 8),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12)
                              ))
    res <- bpi(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_equal(unname(res), c(4, NA, 3))
    spectraData(mse)$basePeakIntensity <- c(5, 0, 4)
    res <- bpi(mse)
    expect_equal(unname(res), c(5, 0, 4))
    res <- bpi(mse, initial = FALSE)
    expect_equal(unname(res), c(4, 0, 3))

    res <- tic(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_equal(unname(res), c(8, 0, 12))
    res <- tic(mse, initial = FALSE)
    expect_equal(unname(res), c(10, 0, 6))
})

test_that("scanIndex,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12)
                              ))
    res <- scanIndex(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_equal(unname(res), c(2L, NA, NA))
    spectraData(mse)$spIdx <- c(3L, 4L, 5L)
    res <- scanIndex(mse)
    expect_equal(unname(res), c(3L, 4L, 5L))
})

test_that("smoothed, smoothed<-,MSnExperiment work", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE)
                              ))
    res <- smoothed(mse)
    expect_equal(names(res), rownames(mse@spectraData))
    expect_equal(unname(res), c(NA, NA, TRUE))
    smoothed(mse) <- TRUE
    expect_equal(unname(smoothed(mse)), c(TRUE, TRUE, TRUE))
    smoothed(mse) <- c(TRUE, FALSE, FALSE)
    expect_equal(unname(smoothed(mse)), c(TRUE, FALSE, FALSE))
    expect_error(smoothed(mse) <- 3)
    expect_error(smoothed(mse) <- c(TRUE, FALSE))
})

test_that("filterAcquisitionNum,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE)
                              ))
    expect_warning(res <- filterAcquisitionNum(mse, n = c(1L, 3L)))
    expect_true(length(res) == 0)
    spectraData(mse)$acquisitionNum <- 1:3
    res <- filterAcquisitionNum(mse, n = c(1L, 3L))
    expect_true(length(res) == 2)
    expect_equal(spectrapply(res), spectrapply(mse)[c(1, 3)])

    expect_error(filterAcquisitionNum(mse, n = "a"))
})

test_that("filterEmptySpectra,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE)
                              ))
    res <- filterEmptySpectra(mse)
    expect_equal(spectrapply(res), spectrapply(mse)[c(1, 3)])
})

test_that("filterMsLevel,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:4, intensity = 1:4,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE)
                              ))
    res <- filterMsLevel(mse)
    expect_equal(spectrapply(res), spectrapply(mse))
    res <- filterMsLevel(mse, msLevel = 1)
    expect_equal(spectrapply(res), spectrapply(mse)[1])
    res <- filterMsLevel(mse, msLevel = 3)
    expect_equal(length(res), 0)
    expect_true(is(res, "MSnExperiment"))
})

test_that("filterMz,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:10, intensity = 1:10,
                                  scanIndex = 2L),
                              new("Spectrum2"),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE)
                              ))
    expect_equal(filterMz(mse), mse)
    res <- filterMz(mse, mz = c(1, 3))
    expect_true(length(res@processingQueue) == 1)
    expect_warning(sps <- as(res, "list"))
    expect_equal(mz(sps[[1]]), 1:3)
    expect_true(isEmpty(sps[[2]]))
    expect_equal(mz(sps[[3]]), 1:3)

    res <- filterMz(mse, mz = c(1, 2), msLevel. = 2)
    expect_warning(sps <- as(res, "list"))
    expect_equal(mz(sps[[1]]), 1:10)
    expect_equal(mz(sps[[3]]), 1:2)

    ## errors
    expect_error(filterMz(mse, mz = 1), "'mz' must be")
    expect_error(filterMz(mse, mz = c("a", "n")), "'mz' must be")
    expect_error(filterMz(mse, mz = c(1, 2), msLevel. = "a"), "'msLevel' must")
})

test_that("filterPrecursorScan,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:10, intensity = 1:10,
                                  scanIndex = 2L, acquisitionNum = 3L),
                              new("Spectrum2", precScanNum = 3L,
                                  acquisitionNum = 4L),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, smoothed = TRUE,
                                  acquisitionNum = 5L, precScanNum = 2L),
                              new("Spectrum2", mz = 1:5, intensity = 1:5,
                                  acquisitionNum = 6L, precScanNum = 3L)
                              ))
    expect_equal(filterPrecursorScan(mse), mse)
    res <- filterPrecursorScan(mse, acquisitionNum = 8L)
    expect_true(length(res) == 0)
    res <- filterPrecursorScan(mse, acquisitionNum = 3L)
    expect_equal(length(res), 3)
    expect_equal(unname(acquisitionNum(res)), c(3L, 4L, 6L))
    res <- filterPrecursorScan(mse, acquisitionNum = 6L)
    expect_equal(unname(acquisitionNum(res)), c(3L, 6L))
    res <- filterPrecursorScan(mse, acquisitionNum = 5L)
    expect_equal(unname(acquisitionNum(res)), 5L)
})

test_that("filterPolarity,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:10, intensity = 1:10,
                                  scanIndex = 2L, acquisitionNum = 3L),
                              new("Spectrum2", precScanNum = 3L,
                                  acquisitionNum = 4L, polarity = 0L),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, polarity = 0L,
                                  acquisitionNum = 5L, precScanNum = 2L),
                              new("Spectrum2", mz = 1:5, intensity = 1:5,
                                  acquisitionNum = 6L, polarity = 1L)
                              ))
    expect_equal(filterPolarity(mse), mse)
    res <- filterPolarity(mse, polarity. = "a")
    expect_equal(length(res), 0)
    res <- filterPolarity(mse, polarity. = 0)
    expect_equal(unname(acquisitionNum(res)), c(4L, 5L))
    res <- filterPolarity(mse, polarity. = 1)
    expect_equal(unname(acquisitionNum(res)), 6L)
})

test_that("filterRt,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:10, intensity = 1:10,
                                  scanIndex = 2L, acquisitionNum = 3L, rt = 1),
                              new("Spectrum2", precScanNum = 3L, rt = 1.1,
                                  acquisitionNum = 4L, polarity = 0L),
                              new("Spectrum2", mz = 1:3, intensity = 1:3,
                                  tic = 12, polarity = 0L, rt = 1.4,
                                  acquisitionNum = 5L, precScanNum = 2L),
                              new("Spectrum2", mz = 1:5, intensity = 1:5,
                                  acquisitionNum = 6L, polarity = 1L, rt = 1.6)
                              ))
    expect_equal(mse, filterRt(mse))
    res <- filterRt(mse, rt = c(2, 4))
    expect_true(length(res) == 0)
    res <- filterRt(mse, rt = c(1, 1.1))
    expect_equal(unname(acquisitionNum(res)), c(3L, 4L))
    res <- filterRt(mse, rt = c(1, 1.1), msLevel. = 2)
    expect_equal(unname(acquisitionNum(res)), c(3L, 4L))
    res <- filterRt(mse, rt = c(1.2, 2))
    expect_equal(unname(acquisitionNum(res)), c(5L, 6L))
    res <- filterRt(mse, rt = c(1.2, 2), msLevel. = 2)
    expect_equal(unname(acquisitionNum(res)), c(3L, 5L, 6L))

    expect_error(filterRt(mse, rt = 4), "'rt' must be")
    expect_error(filterRt(mse, rt = 1:4), "'rt' must be")
    expect_error(filterRt(mse, rt = c(FALSE, TRUE)), "'rt' must be")
    expect_error(filterRt(mse, rt = c(1, 2), msLevel. = "z"), "'msLevel'")
})

test_that("$,MSnExperiment works", {
    spl <- list(new("Spectrum1", mz = 1:10, intensity = 1:10,
                    fromFile = 1L, scanIndex = 1L),
                new("Spectrum2", fromFile = 2L, scanIndex = 1L),
                new("Spectrum2", mz = 1:3, intensity = 1:3,
                    fromFile = 1L, scanIndex = 3L),
                new("Spectrum2", mz = 1:5, intensity = 1:5,
                    fromFile = 2L, scanIndex = 4L)
                )
    mse <- MSnExperiment(spl, sampleData = DataFrame(sample_idx = 1:2,
                                                     sample_name = letters[1:2]))
    expect_equal(colnames(sampleData(mse)), c("sample_idx", "sample_name"))
    expect_equal(mse$sample_idx, 1:2)
    expect_equal(mse$sample_name, c("a", "b"))
    mse$sample_name <- c("c", "z")
    expect_equal(mse$sample_name, c("c", "z"))
    mse$sample_desc <- c("c", "z")
    expect_equal(ncol(sampleData(mse)), 3)

    expect_error(mse$sample_name <- c("c", "z", "b"))
})

test_that("splitByFile,MSnExperiment works", {
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
           system.file("microtofq/MM8.mzML", package = "msdata"))
    inMem <- readMSnExperiment(f, backend = BackendMemory())
    expect_error(splitByFile(inMem, f = factor(1:3)))
    spl <- splitByFile(inMem, f = factor(c("b", "a")))
    expect_equal(spectrapply(spl[[1]]), spectrapply(filterFile(inMem, 2)))
    expect_equal(sampleData(spl[[1]]), sampleData(filterFile(inMem, 2)))
    expect_equal(spectrapply(spl[[2]]), spectrapply(filterFile(inMem, 1)))
    expect_equal(sampleData(spl[[2]]), sampleData(filterFile(inMem, 1)))
})

test_that("bin,MSnExperiment works", {
    mse <- MSnExperiment(list(new("Spectrum1", mz = 1:10, intensity = 1:10),
                              new("Spectrum2", mz = 1:10, intensity = 11:20)
                              ))
    expect_warning(expect_equal(mse, bin(mse, msLevel. = 3)))
    res <- bin(mse, binSize = 2)
    sps <- as(res, "list")
    expect_equal(mz(sps[[1]]), c(2, 4, 6, 8, 10))
    expect_equal(intensity(sps[[1]]), c(3, 7, 11, 15, 19))
    expect_equal(mz(sps[[2]]), c(2, 4, 6, 8, 10))

    res <- bin(mse, binSize = 2, msLevel. = 2)
    sps <- as(res, "list")
    expect_equal(mz(sps[[1]]), 1:10)
    expect_equal(intensity(sps[[1]]), 1:10)
    expect_equal(mz(sps[[2]]), c(2, 4, 6, 8, 10))
})

test_that("estimateNoise,MSnExperiment works", {
    res <- estimateNoise(sciex_inmem[1:10], method = "SuperSmoother")
    expect_equal(res[[1]],
                 estimateNoise(sciex_inmem[[1]], method = "SuperSmoother"))
})

test_that("compareSpectra,MSnExperiment works", {
    a <- sciex[1:10]
    b <- sciex_inmem[1:10]
    res_a <- compareSpectra(a)
    res_b <- compareSpectra(b)
    expect_equal(res_a, res_b)
    expect_equal(dim(res_b), c(10, 10))
    expect_equal(rownames(res_b), featureNames(b))
    expect_equal(res_b[1, 2], compareSpectra(b[[1]], b[[2]]))

    res_b <- compareSpectra(b, fun = "dotproduct")
    expect_equal(res_b[1, 4], compareSpectra(b[[1]], b[[4]], fun = "dotproduct"))
})

test_that("estimateMzResolution,MSnExperiment works", {
    res <- estimateMzResolution(sciex_inmem[1:10])
    expect_equal(length(res), 10)
    expect_equal(names(res), featureNames(sciex_inmem)[1:10])
    expect_equal(res[[1]], estimateMzResolution(sciex_inmem[[1]]))
})

test_that("normalize,MSnExperiment works", {
    res <- normalize(sciex_inmem[1:10])
    expect_equal(length(res@processingQueue), 1)
    expect_equal(res[[1]], normalize(sciex_inmem[[1]]))
})

test_that("pickPeaks,MSnExperiment works", {
    res <- pickPeaks(sciex_inmem[1:10], method = "SuperSmoother",
                     refineMz = "descendPeak", signalPercentage = 45)
    expect_equal(length(res@processingQueue), 1)
    sps <- as(res, "list")
    expect_equal(sps[[2]], pickPeaks(sciex_inmem[[2]], method = "SuperSmoother",
                                     refineMz = "descendPeak",
                                     signalPercentage = 45))
})

## test_that("quantify,MSnExperiment works", {
##     x1 <- extdata_mzXML_in_mem_ms2
##     f <- dir(system.file(package = "MSnbase", dir = "extdata"),
##              full.name = TRUE, pattern = "mzXML$")
##     x2 <- readMSnExperiment(f)
##     e1 <- quantify(x1, method = "max", reporters = iTRAQ4, verbose = FALSE)
##     e2 <- quantify(x2, method = "max", reporters = iTRAQ4, verbose = FALSE)
##     expect_identical(exprs(e1), exprs(e2))

## })

test_that("removeReporters,MSnExperiment works", {
    in_mem <- tmt_erwinia_in_mem_ms2
    in_mem_rem <- removeReporters(in_mem, TMT6)
    f <- msdata::proteomics(
                     full.names = TRUE,
                     pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    od <- readMSnExperiment(f)
    res <- removeReporters(filterMsLevel(od, 2L), TMT6)
    expect_true(length(res@processingQueue) == 1)

    expect_identical(unname(spectra(in_mem_rem)), unname(as(res, "list")))

    ## Do the call on the full data set.
    res <- removeReporters(od, TMT6)
    expect_identical(unname(spectra(in_mem_rem)),
                     unname(as(filterMsLevel(res, 2L), "list")))

    expect_error(removeReporters(sciex_mzr, TMT6), "No MS level > 1")
})

test_that("smooth,MSnExperiment works", {
    res <- smooth(sciex_mzr[1:10], method = "MovingAverage", halfWindowSize = 3L)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(as(res, "list")[[3]],
                 smooth(sciex_mzr[[3]], method = "MovingAverage",
                        halfWindowSize = 3L))
})

test_that("chromatogram,MSnExperiment works", {
    f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
           system.file("microtofq/MM8.mzML", package = "msdata"))
    mem <- readMSnExperiment(f, backend = BackendMemory())
    mzr <- readMSnExperiment(f, backend = BackendMzR())
    ## file 1: rt 270-307, mz = 94.8, 1004
    ## file 2: rt 0.4-66.7 mz = 95, 1005
    res <- chromatogram(mem, mz = c(1, 2), rt = c(400, 402))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == 2)

    mzm <- matrix(c(100, 120, 200, 220, 300, 320), nrow = 3, byrow = TRUE)
    rtm <- matrix(c(50, 300), nrow = 1)
    res <- chromatogram(mem, mz = mzm, rt = rtm)

    ## Check that the values for all ranges is within the specified ranges
    for (i in 1:nrow(mzm)) {
        expect_true(all(mz(res[i, 1]) >= mzm[i, 1] &
                        mz(res[i, 1]) <= mzm[i, 2]))
        expect_true(all(mz(res[i, 2]) >= mzm[i, 1] &
                        mz(res[i, 2]) <= mzm[i, 2]))
        expect_true(all(rtime(res[i, 1]) >= rtm[1, 1] &
                        rtime(res[i, 1]) <= rtm[1, 2]))
        expect_true(all(rtime(res[i, 2]) >= rtm[1, 1] &
                        rtime(res[i, 2]) <= rtm[1, 2]))
    }
    ## Check that values are correct.
    flt <- filterMz(filterRt(mem, rt = rtm[1, ]), mz = mzm[2, ])
    ints <- split(unlist(lapply(as(flt, "list"), function(z) sum(intensity(z)))),
                  fromFile(flt))
    expect_equal(ints[[1]], intensity(res[2, 1]))
    expect_equal(ints[[2]], intensity(res[2, 2]))
    expect_equal(split(rtime(flt), fromFile(flt))[[1]], rtime(res[2, 1]))
    expect_equal(split(rtime(flt), fromFile(flt))[[2]], rtime(res[2, 2]))
    ## fData
    expect_true(nrow(fData(res)) == nrow(res))
    expect_true(all(colnames(fData(res)) == c("mzmin", "mzmax",
                                              "rtmin", "rtmax", "polarity")))
    expect_true(all(fData(res)$rtmin == 50))
    expect_true(all(fData(res)$rtmax == 300))
    expect_equal(fData(res)$mzmin, c(100, 200, 300))
    expect_equal(fData(res)$mzmax, c(120, 220, 320))
    expect_equal(fData(res)$polarity, c(1, 1, 1))

    ## Now with ranges for which we don't have values in one or the other.
    rtr <- matrix(c(280, 300, 20, 40), nrow = 2,
                  byrow = TRUE)  ## Only present in first, or 2nd file
    res <- chromatogram(mzr, rt = rtr)
    ## Check fromFile
    for (i in 1:ncol(res))
        expect_true(all(sapply(res[, i], fromFile) == i))
    expect_equal(length(res[2, 1]), 0)
    expect_equal(length(res[1, 2]), 0)
    ## Check rtime
    expect_true(all(rtime(res[1, 1]) >= rtr[1, 1] &
                    rtime(res[1, 1]) <= rtr[1, 2]))
    expect_true(all(rtime(res[2, 2]) >= rtr[2, 1] &
                    rtime(res[2, 2]) <= rtr[2, 2]))
    ## Check intensity
    flt <- filterRt(mzr, rt = rtr[1, ])
    spctr <- split(as(flt, "list"), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[1, 1]))
    flt <- filterRt(mzr, rt = rtr[2, ])
    spctr <- split(as(flt, "list"), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[2, 2]))

    ## Check that phenoType is correctly passed.
    pd <- data.frame(name = c("first", "second"), idx = 1:2)
    sampleData(mzr) <- DataFrame(pd)
    chrs <- chromatogram(mzr)
    ## rownames(pd) <- colnames(chrs)
    expect_equal(pData(chrs), pd)

    chrs_2 <- chromatogram(mzr, msLevel = 1:4)
    expect_equal(chrs, chrs_2)
})
