context("MSnExp filter functions")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem2 <- readMSData(f, verbose = FALSE)  ## That's the MS 2 data.
inmem1 <- readMSData(f, verbose = FALSE, msLevel = 1)  ## MS 1 data.
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
