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
    tmp <- removePeaks(sciex_inmem[900:920], t = 5000)
    tmp <- applyProcessingQueue(tmp)
    expect_equal(tmp@processingQueue, list())
    expect_equal(tmp@backend@modCount, 1L)
    expect_error(setBackend(tmp, BackendMzR()), "Can not change backend")

    ## Memory, modify -> Hdf5
    tmp_h5 <- setBackend(tmp, BackendHdf5(),
                         path = paste0(tempdir(), "/switch2/"))
    expect_true(length(tmp_h5@processingQueue) == 0)
    expect_equal(unname(tmp_h5@backend@modCount), 1L)
    sps_h5 <- spectrapply(tmp_h5)
    sps <- spectrapply(sciex_inmem[900:920])
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
