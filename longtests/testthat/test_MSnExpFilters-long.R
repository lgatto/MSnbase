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
