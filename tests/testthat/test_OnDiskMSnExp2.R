context("OnDiskMSnExp class")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

test_that("OnDiskMSnExp constructor", {
    x <- readMSData2(f)
    expect_true(validObject(x))
    expect_true(all(unique(msLevel(x)) == 1:2))
    expect_true(all(isCurrent(x)))
    expect_true(isVersioned(x))
    expect_null(show(x))
    x1 <- readMSData2(f, msLevel = 1)
    expect_true(validObject(x1))
    expect_true(unique(msLevel(x1)) == 1)
    x2 <- readMSData2(f, msLevel = 2)
    expect_true(validObject(x2))
    expect_true(unique(msLevel(x2)) == 2)
    expect_identical(length(x), length(x1) + length(x2))
    expect_identical(featureNames(x),
                     sort(c(featureNames(x1), featureNames(x2))))
    expect_identical(fileNames(x), fileNames(x1))
    expect_identical(fileNames(x), fileNames(x2))
})

test_that("compare MS2 on disk and in memory", {
    x1 <- readMSData(f, verbose = FALSE, centroided. = FALSE)
    x2 <- readMSData2(f, msLevel. = 2, centroided. = FALSE)
    expect_identical(length(x1), length(x2))
    expect_false(identical(featureNames(x1), featureNames(x2)))
    featureNames(x2) <- featureNames(x1)
    expect_identical(featureNames(x1), featureNames(x2))
    ## testing all accessors
    expect_identical(x1[[1]], x2[[1]])
    i <- sample(length(x1), 1)
    expect_identical(x1[[i]], x2[[i]])
    ## Test [ with an all.equal method - see issue #122
    expect_identical(abstract(x1), abstract(x2))
    expect_identical(acquisitionNum(x1), acquisitionNum(x2))
    expect_identical(analyser(x1), analyser(x2))
    expect_identical(analyserDetails(x1), analyserDetails(x2))
    expect_identical(analyzer(x1), analyzer(x2))
    expect_identical(analyzerDetails(x1), analyzerDetails(x2))
    expect_identical(collisionEnergy(x1), collisionEnergy(x2))
    expect_identical(description(x1), description(x2))
    expect_identical(detectorType(x1), detectorType(x2))
    expect_identical(dim(x1), dim(x2))
    expect_identical(estimateNoise(x1[[1]]), estimateNoise(x2[[1]]))
    expect_identical(expemail(x1), expemail(x2))
    expect_identical(experimentData(x1), experimentData(x2))
    expect_identical(expinfo(x1), expinfo(x2))
    expect_identical(exptitle(x1), exptitle(x2))
    expect_identical(fileNames(x1), fileNames(x2))
    expect_identical(fromFile(x1), fromFile(x2))
    expect_identical(instrumentCustomisations(x1), instrumentCustomisations(x2))
    expect_identical(instrumentManufacturer(x1), instrumentManufacturer(x2))
    expect_identical(instrumentModel(x1), instrumentModel(x2))
    expect_identical(intensity(x1), intensity(x2))
    expect_identical(ionCount(x1), ionCount(x2))
    expect_identical(ionSource(x1), ionSource(x2))
    expect_identical(ionSourceDetails(x1), ionSourceDetails(x2))
    expect_identical(msLevel(x1), msLevel(x2))
    expect_identical(notes(x1), notes(x2))
    expect_identical(pData(x1), pData(x2))
    expect_identical(phenoData(x1), phenoData(x2))
    expect_identical(polarity(x1), polarity(x2))
    expect_identical(precScanNum(x1), precScanNum(x2))
    expect_identical(precursorIntensity(x1), precursorIntensity(x2))
    expect_identical(protocolData(x1), protocolData(x2))
    expect_identical(pubMedIds(x1), pubMedIds(x2))
    pubMedIds(x1) <- pubMedIds(x2) <- "23692960"
    expect_identical(pubMedIds(x1), pubMedIds(x2))
    expect_identical(rtime(x1), rtime(x2))
    expect_identical(sampleNames(x1), sampleNames(x2))
    sampleNames(x1) <- sampleNames(x2) <- "X1"
    expect_identical(sampleNames(x1), sampleNames(x2))
    expect_identical(scanIndex(x1), scanIndex(x2))
    expect_identical(spectra(x1), spectra(x2))
    expect_identical(tic(x1), tic(x2))
})

test_that("Default and setting centroided", {
    x1 <- readMSData(f, verbose = FALSE, centroided. = NA)
    x2 <- readMSData2(f, msLevel. = 2, centroided. = NA)
    featureNames(x2) <- featureNames(x1)
    expect_identical(centroided(x1), centroided(x2))
    x1 <- readMSData(f, verbose = FALSE, centroided. = TRUE)
    x2 <- readMSData2(f, msLevel. = 2, centroided. = TRUE)
    featureNames(x2) <- featureNames(x1)
    expect_identical(centroided(x1), centroided(x2))
    centroided(x2) <- centroided(x1) <- FALSE
    expect_identical(centroided(x1), centroided(x2))
})

test_that("Write mgf", {
    x1 <- readMSData(f, verbose = FALSE)
    x2 <- readMSData2(f, msLevel. = 2)
    tf1 <- tempfile()
    tf2 <- tempfile()
    writeMgfData(x1[2:4], con = tf1)
    writeMgfData(x2[2:4], con = tf2)
    mgf1 <- readMgfData(tf1)
    mgf2 <- readMgfData(tf2)
    mgf2@processingData <- mgf1@processingData
    mgf2@phenoData <- mgf1@phenoData
    expect_equal(mgf1, mgf2)
    unlink(tf1)
    unlink(tf2)
})

test_that("Adding identification data", {
    quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
    x1 <- readMSData(quantFile, verbose = FALSE)
    x1 <- addIdentificationData(x1, identFile, verbose = FALSE)
    x2 <- readMSData2(quantFile, verbose = FALSE)
    x2 <- addIdentificationData(x2, identFile, verbose = FALSE)
    fv <- intersect(fvarLabels(x1), fvarLabels(x2))
    expect_identical(fData(x1)[, fv], fData(x2)[, fv])
    expect_identical(idSummary(x1), idSummary(x2))
    x1rm <- removeMultipleAssignment(x1)
    x2rm <- removeMultipleAssignment(x2)
    expect_identical(fData(x1rm)[, fv], fData(x2rm)[, fv])
    x1nid <- removeNoId(x1)
    x2nid <- removeNoId(x2)
    expect_identical(fData(x1nid)[, fv], fData(x2nid)[, fv])
})


test_that("Spectrum processing", {
    ## bin
    ## clean -> test_OnDiskMSnExp.R
    ## normalise/normalize -> test_OnDiskMSnExp_other_methods.R
    ## pickPeaks -> test_OnDiskMSnExp_other_methods.R
    ## quantify
    ## removePeaks -> test_OnDiskMSnExp.R
    ## removeReporters -> here
    ## smooth -> test_OnDiskMSnExp_other_methods.R
    ## trimMz -> test_OnDiskMSnExp_other_methods.R
})

test_that("removeReporters,OnDiskMSnExp", {
    in_mem <- readMSData(f, msLevel. = 2)
    in_mem_rem <- removeReporters(in_mem, TMT6)

    on_disk <- readMSData2(f)
    on_disk_2 <- filterMsLevel(on_disk, msLevel. = 2)
    on_disk_2_rem <- removeReporters(on_disk_2, TMT6)
    sp_rem <- spectra(on_disk_2_rem)
    expect_identical(unname(spectra(in_mem_rem)), unname(sp_rem))

    ## Do the call on the full data set.
    on_disk_rem <- removeReporters(on_disk, TMT6)
    sp_rem <- spectra(on_disk_rem)
    ## Subset to MS2 spectra
    sp_rem <- sp_rem[unlist(lapply(sp_rem, msLevel)) > 1]
    expect_identical(unname(spectra(in_mem_rem)), unname(sp_rem))
})
