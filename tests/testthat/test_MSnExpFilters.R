context("MSnExp filter functions")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem <- readMSData(f, verbose = FALSE)  ## That's the MS 2 data.
inmem1 <- readMSData(f, verbose = FALSE, msLevel = 1)  ## MS 1 data.
ondisk <- MSnbase:::readMSData2(f, verbose = FALSE)
ondisk2 <- MSnbase:::readMSData2(f, msLevel = 2, verbose = FALSE)

test_that("filterMsLevel", {
    expect_true(all.equal(inmem, filterMsLevel(inmem, msLevel. = 2)))
    expect_true(all.equal(ondisk2, filterMsLevel(ondisk2, msLevel. = 2)))
    ## Compare in-mem MS2 with filterMsLevel of on-disk MS 1 and 2.
    expect_true(all.equal(inmem, filterMsLevel(ondisk, msLevel. = 2)))
    ## Compare in-mem MS1 with filterMsLevel of on-disk MS 1 and 2
    expect_true(all.equal(inmem1, filterMsLevel(ondisk, msLevel. = 1)))
    ## in-mem MS 2 with on-disk MS 2
    expect_true(all.equal(inmem, ondisk2))
    expect_true(all.equal(filterMsLevel(inmem, msLevel. = 2),
                          filterMsLevel(ondisk2, msLevel. = 2)))
    expect_identical(length(filterMsLevel(inmem, 1)), 0L)
    expect_identical(length(filterMsLevel(ondisk2, 1)), 0L)
    expect_true(all.equal(filterMsLevel(ondisk, 2), ondisk2))
    expect_identical(length(filterMsLevel(ondisk, 1)),
                     sum(fData(ondisk)$msLevel == 1))
})
