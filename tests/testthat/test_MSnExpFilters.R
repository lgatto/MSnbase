context("MSnExp filter functions")

inmem1 <- tmt_erwinia_in_mem_ms1
inmem2 <- tmt_erwinia_in_mem_ms2
ondisk <- tmt_erwinia_on_disk
ondisk1 <- tmt_erwinia_on_disk_ms1
ondisk2 <- tmt_erwinia_on_disk_ms2


test_that("filterRt", {
    rtlim <- c(122, 1300)
    expect_true(all.equal(filterRt(inmem1, rtlim, 1),
                          filterRt(ondisk1, rtlim, 1)))
    expect_true(all.equal(filterRt(inmem2, rtlim, 2),
                          filterRt(ondisk2, rtlim, 2)))
    expect_true(all.equal(filterRt(inmem2, rtlim, 2),
                          filterMsLevel(filterRt(ondisk, rtlim, 2), 2)))
})

test_that("filterFile", {
    ## Use two files.
    mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                 system.file("microtofq/MM8.mzML", package = "msdata"))
    oneFileInMem <- readMSData(mzfiles[2], verbose = FALSE, msLevel = 1)
    twoFileInMem <- microtofq_in_mem_ms1
    twoFileOnDisk <- microtofq_on_disk
    secondFileOnDisk <- readMSData(mzfiles[2], verbose = FALSE, mode = "onDisk")
    ## Note: all.equal MSnExp, MSnExp will fail because of the
    ## experimentData and featureNames
    expect_true(all.equal(spectra(filterFile(twoFileInMem, file = 2)),
                          spectra(oneFileInMem), check.names = FALSE))
    expect_true(all.equal(filterFile(twoFileOnDisk, file = 2),
                          oneFileInMem, check.names = FALSE))
    ## Below breaks because of fromFile
    expect_true(all.equal(filterFile(twoFileOnDisk, file = 2),
                          filterFile(twoFileInMem, file = 2),
                          check.names = FALSE))
    ## Check experimentData:
    secondFile <- filterFile(twoFileOnDisk, file = 2)
    expect_equal(experimentData(secondFile)@instrumentManufacturer,
                 experimentData(oneFileInMem)@instrumentManufacturer)
    ## Compare the sub set to the single file.
    expect_true(all.equal(secondFileOnDisk, filterFile(twoFileOnDisk, file = 2)))
    expect_identical(pData(secondFileOnDisk),
                     pData(filterFile(twoFileOnDisk, file = 2)))
    expect_identical(experimentData(secondFileOnDisk),
                     experimentData(filterFile(twoFileOnDisk, file = 2)))
})

test_that("filterAcquisitionNum", {
    filtered <- filterAcquisitionNum(ondisk, n = 1000:1100)
    expect_true(all(acquisitionNum(filtered) %in% 1000:1100))
    spctr <- spectra(filtered)
    expect_true(all(unlist(lapply(spctr, acquisitionNum)) %in% 1000:1100))
    ## ## Use unavailable acquisition numbers
    ## expect_warning(ll <- length(filterAcquisitionNum(ondisk, n = 1:100)))
    ## expect_true(ll == 0)
    ## Compare on-disk with in-mem.
    expect_true(all.equal(filterAcquisitionNum(ondisk1, n = 1000:1100),
                          filterAcquisitionNum(inmem1, n = 1000:1100)))
    expect_true(all.equal(filterAcquisitionNum(ondisk2, n = 1000:1100),
                          filterAcquisitionNum(inmem2, n = 1000:1100)))
    ## Torture tests. The two files have different number of spectra.
    mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                 system.file("microtofq/MM8.mzML", package = "msdata"))
    twoFileOnDisk <- microtofq_on_disk
    centroided(twoFileOnDisk) <- TRUE
    secondFile <- readMSData(mzfiles[2], verbose = FALSE, centroided. = TRUE, mode = "onDisk")
    expect_warning(res <- filterAcquisitionNum(twoFileOnDisk, n = 180:190, file = 1))
    expect_identical(fileNames(res), fileNames(twoFileOnDisk)[2])
    ## contains basically only the second file.
    expect_true(all.equal(secondFile, res))
    expect_identical(pData(secondFile), pData(res))
    expect_identical(experimentData(secondFile), experimentData(res))
    ## Second subset.
    res <- filterAcquisitionNum(twoFileOnDisk, n = 180:190)
    expect_true(all(acquisitionNum(res) %in% 180:190))
    expect_identical(pData(secondFile), pData(res))
    expect_identical(experimentData(secondFile), experimentData(res))
})


test_that("filterPrecursorScan", {
    expect_error(filterPrecursorScan(inmem1, 1),
                 "column\\(s\\) acquisitionNum/precursorScanNum is/are missing")
    expect_true(all.equal(ondisk[1003:1013], filterPrecursorScan(ondisk, 1003)))
    expect_true(all.equal(ondisk[1], filterPrecursorScan(ondisk, 1)))
})


test_that("filterPolarity", {
    fls <- dir(system.file("sciex", package = "msdata"),
               full.names = TRUE, recursive = TRUE)
    rw <- readMSData(files = fls[1], mode = "onDisk")
    expect_identical(length(filterPolarity(rw, -1)), 0L)
    expect_identical(length(filterPolarity(rw, 1)), length(rw))
    pol2 <- rep(c(1, -1), length.out = length(rw))
    fData(rw)$polarity <- pol2
    expect_identical(length(filterPolarity(rw, 1)), sum(pol2 == 1))
    expect_identical(length(filterPolarity(rw, -1)), sum(pol2 == -1))
})

test_that("filterPrecursorMz works", {
    res <- filterPrecursorMz(tmt_im_ms2_sub, mz = 417.75, ppm = 20)
    expect_true(length(res) == 1)

    res <- filterPrecursorMz(tmt_im_ms2_sub, mz = 417.75, ppm = 100)
    expect_true(length(res) == 2)

    res <- filterPrecursorMz(tmt_od_sub, mz = 417.75, ppm = 20)
    expect_true(length(res) == 1)
    res <- filterPrecursorMz(tmt_od_sub, mz = 417.75, ppm = 100)
    expect_true(length(res) == 2)

    res <- filterPrecursorMz(tmt_od_sub)
    expect_equal(length(res), length(tmt_od_sub))

    expect_error(filterPrecursorMz(sciex, mz = 233), "MS1 spectra")

    expect_error(filterPrecursorMz(tmt_od_sub, mz = c(1, 2)), "single m/z")
})

test_that("filterIsolationWindow works", {
    res <- filterIsolationWindow(tmt_od_sub)
    expect_equal(length(res), length(tmt_od_sub))

    res <- filterIsolationWindow(tmt_od_sub, mz = 411)
    expect_equal(length(res), 2)
    expect_true(all(isolationWindowLowerMz(res) < 411))
    expect_true(all(isolationWindowUpperMz(res) > 411))

    expect_error(filterIsolationWindow(tmt_od_sub, mz = c(1, 2)), "single m/z")

    expect_error(filterIsolationWindow(tmt_im_ms2_sub, mz = 411), "not available")

    res <- filterIsolationWindow(sciex, mz = 411)
    expect_true(length(res) == 0)
})
