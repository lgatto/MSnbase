context("MSnExp filter functions")

inmem1 <- tmt_erwinia_in_mem_ms1
inmem2 <- tmt_erwinia_in_mem_ms2
ondisk <- tmt_erwinia_on_disk
ondisk1 <- tmt_erwinia_on_disk_ms1
ondisk2 <- tmt_erwinia_on_disk_ms2

test_that("filterMsLevel", {
    expect_true(all.equal(inmem2, filterMsLevel(inmem2, msLevel. = 2)))
    expect_true(all.equal(ondisk2, filterMsLevel(ondisk2, msLevel. = 2)))
    ## Compare in-mem MS2 with filterMsLevel of on-disk MS 1 and 2.
    expect_true(all.equal(inmem2, filterMsLevel(ondisk, msLevel. = 2)))
    ## Compare in-mem MS1 with filterMsLevel of on-disk MS 1 and 2
    expect_true(all.equal(inmem1, filterMsLevel(ondisk, msLevel. = 1)))
    ## in-mem MS 2 with on-disk MS 2
    expect_true(all.equal(inmem2, ondisk2))
    expect_true(all.equal(filterMsLevel(inmem2, msLevel. = 2),
                          filterMsLevel(ondisk2, msLevel. = 2)))
    expect_identical(length(filterMsLevel(inmem2, 1)), 0L)
    expect_identical(length(filterMsLevel(ondisk2, 1)), 0L)
    expect_true(all.equal(filterMsLevel(ondisk, 2), ondisk2))
    expect_identical(length(filterMsLevel(ondisk, 1)),
                     sum(fData(ondisk)$msLevel == 1))
})

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
    ## Use unavailable acquisition numbers
    expect_warning(ll <- length(filterAcquisitionNum(ondisk, n = 1:100)))
    expect_true(ll == 0)
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
    secondFile <- readMSData(mzfiles[2], verbose = FALSE, centroided = TRUE, mode = "onDisk")
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

test_that("filterMz", {
    ## Subset all objects to speed up tests.
    ondisk1 <- filterRt(ondisk1, rt = c(1150, 1200))
    ondisk2 <- filterRt(ondisk2, rt = c(1150, 1200))
    ondisk <- filterRt(ondisk, rt = c(1150, 1200))
    inmem1 <- filterRt(inmem1, rt = c(1150, 1200))
    inmem2 <- filterRt(inmem2, rt = c(1150, 1200))
    ## On single, multi MS file.
    ## o Test for MSnExp, MS 1 and 2
    inMemF <- filterMz(inmem1, mz = c(500, 550))
    mzr <- range(mz(inMemF))
    expect_true(mzr[1] >= 500 & mzr[2] <= 550)
    inMemF <- filterMz(inmem2, mz = c(300, 900))
    mzr <- range(mz(inMemF))
    expect_true(mzr[1] >= 300 & mzr[2] <= 900)
    ## o Test for OnDiskMSnExp, MS 1, MS 2 and both.
    onDisk1F <- filterMz(ondisk1, mz = c(500, 550))
    mzr <- range(mz(onDisk1F))
    expect_true(mzr[1] >= 500 & mzr[2] <= 550)
    onDisk2F <- filterMz(ondisk2, mz = c(300, 900))
    mzr <- range(mz(onDisk2F))
    expect_true(mzr[1] >= 300 & mzr[2] <= 900)
    onDiskF_all <- filterMz(ondisk, mz = c(500, 550))
    expect_true(all.equal(onDisk1F, filterMsLevel(onDiskF_all, msLevel. = 1)))
    ##  Do the filtering separately for MS1 and MS2 levels:
    onDiskF <- filterMz(ondisk, mz = c(500, 550), msLevel. = 1)
    onDiskF <- filterMz(onDiskF, mz = c(300, 900), msLevel. = 2)
    expect_true(all.equal(onDisk1F, filterMsLevel(onDiskF, msLevel. = 1)))
    expect_true(all.equal(onDisk2F, filterMsLevel(onDiskF, msLevel. = 2)))
    ##  What if we're asking for an msLevel that doesn't exist.
    Test <- filterMz(ondisk1, mz = c(500, 550), msLevel. = 2)
    expect_true(all.equal(Test, ondisk1)) ## no change
    expect_false(isTRUE(all.equal(Test, onDisk1F)))
    expect_is(all.equal(Test, onDisk1F), "character")
    ## o Compare between MSnExp and OnDiskMSnExp.
    onDiskF <- filterMz(ondisk1, mz = c(500, 550))
    inMemF <- filterMz(inmem1, mz = c(500, 550))
    expect_true(all.equal(onDiskF, inMemF))
    onDiskF <- filterMz(ondisk2, mz = c(300, 900))
    inMemF <- filterMz(inmem2, mz = c(300, 900))
    expect_true(all.equal(onDiskF, inMemF))
    ## On multiple files.
    mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                 system.file("microtofq/MM8.mzML", package = "msdata"))
    twoFileOnDisk <- microtofq_on_disk
    twoFileOnDiskF <- filterMz(twoFileOnDisk, mz = c(300, 350))
    mzr <- range(mz(twoFileOnDiskF))
    expect_true(mzr[1] >= 300 & mzr[2] <= 350)
})

test_that("filterPrecursorScan", {
    expect_error(filterPrecursorScan(inmem1, 1),
                 "column\\(s\\) acquisitionNum/precursorScanNum is/are missing")
    expect_true(all.equal(ondisk[1:11], filterPrecursorScan(ondisk, 1003)))
    expect_true(all.equal(ondisk[c(1:11, 19:20)],
                          filterPrecursorScan(ondisk, c(1003, 1022))))
    expect_true(all.equal(ondisk[NA], filterPrecursorScan(ondisk, 1)))
})
