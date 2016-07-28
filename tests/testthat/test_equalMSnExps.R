context("[OnDisk]MSnExp equality")

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
## Set centroided to NA to fit with default behaviour of readMSData2
inmem1 <- readMSData(f, msLevel = 1, verbose = FALSE, centroided = NA)
inmem2 <- readMSData(f, msLevel = 2, verbose = FALSE, centroided = NA)
ondisk <- readMSData2(f, verbose = FALSE)
ondisk1 <- readMSData2(f, msLevel = 1, verbose = FALSE)
ondisk2 <- readMSData2(f, msLevel = 2, verbose = FALSE)

test_that("Equality function", {
    ## testing it on spectra only.
    expect_true(all.equal(inmem1, ondisk1))
    expect_true(all.equal(inmem2, ondisk2))
    ## postive controls
    expect_true(all.equal(inmem1, inmem1))
    expect_true(all.equal(inmem2, inmem2))
    expect_true(all.equal(ondisk, ondisk))
    ## negative controls
    expect_false(isTRUE(all.equal(ondisk, ondisk1)))
    expect_false(isTRUE(all.equal(ondisk, ondisk2)))
})
