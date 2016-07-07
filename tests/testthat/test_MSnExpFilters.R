context("MSnExp filter functions")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem2 <- readMSData(f, centroided = NA, verbose = FALSE)  ## That's the MS 2 data.
inmem1 <- readMSData(f, centroided = NA, verbose = FALSE, msLevel = 1)  ## MS 1 data.
ondisk <- readMSData2(f, verbose = FALSE)
ondisk1 <- readMSData2(f, msLevel = 1, verbose = FALSE)
ondisk2 <- readMSData2(f, msLevel = 2, verbose = FALSE)

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
    oneFileInMem <- readMSData(mzfiles[2], verbose = FALSE, msLevel = 1,
                               centroided = TRUE)
    suppressWarnings(
        twoFileInMem <- readMSData(mzfiles, verbose = FALSE, msLevel = 1,
                                   centroided = TRUE)
    )
    suppressWarnings(
        twoFileOnDisk <- readMSData2(mzfiles, verbose = FALSE,
                                     centroided = TRUE)
    )
    secondFileOnDisk <- readMSData2(mzfiles[2], verbose = FALSE,
                                    centroided = TRUE)
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
    suppressWarnings(
        expect_true(length(filterAcquisitionNum(ondisk, n = 1:100)) == 0)
    )
    ## Compare on-disk with in-mem.
    expect_true(all.equal(filterAcquisitionNum(ondisk1, n = 1000:1100),
                          filterAcquisitionNum(inmem1, n = 1000:1100)))
    expect_true(all.equal(filterAcquisitionNum(ondisk2, n = 1000:1100),
                          filterAcquisitionNum(inmem2, n = 1000:1100)))

    ## Torture tests. The two files have different number of spectra.
    mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                 system.file("microtofq/MM8.mzML", package = "msdata"))
    suppressWarnings(
        twoFileOnDisk <- readMSData2(mzfiles, verbose = FALSE,
                                     centroided = TRUE)
    )
    secondFile <- readMSData2(mzfiles[2], verbose = FALSE, centroided = TRUE)
    suppressWarnings(
        res <- filterAcquisitionNum(twoFileOnDisk, n = 180:190, file = 1)
    )
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
