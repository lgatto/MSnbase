context("MSnExp filter functions")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem <- readMSData(f, verbose = FALSE)
ondisk <- MSnbase:::readMSData2(f)
ondisk2 <- MSnbase:::readMSData2(f, msLevel = 2)

test_that("filterMsLevel", {
    expect_true(all.equal(inmem, filterMsLevel(inmem, msLevel. = 2)))
    expect_true(all.equal(ondisk2, filterMsLevel(ondisk2, msLevel. = 2)))
    expect_true(all.equal(inmem, ondisk2))
    expect_true(all.equal(filterMsLevel(inmem, msLevel. = 2),
                          filterMsLevel(ondisk2, msLevel. = 2)))
    expect_identical(length(filterMsLevel(inmem, 1)), 0L)
    expect_identical(length(filterMsLevel(ondisk2, 1)), 0L)
    expect_true(all.equal(filterMsLevel(ondisk, 2), ondisk2))
    expect_identical(length(filterMsLevel(ondisk, 1)),
                     sum(fData(ondisk)$msLevel == 1))
})
