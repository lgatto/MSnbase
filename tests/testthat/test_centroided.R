test_that("centroided accessor with/without na.fail", {
    x <- tmt_erwinia_in_mem_ms2
    x2 <- tmt_erwinia_on_disk
    ## in-mem, single spectrum
    expect_true(centroided(x[[1]]))
    ## in-mem, experiment
    expect_true(all(centroided(x)))
    ## on-disk, single spectrum
    expect_true(!centroided(x2[[1]]))
    expect_true(centroided(x2[[2]]))
    ## on-disk, experiment
    expect_true(all(!centroided(filterMsLevel(x2, 1))))
    expect_true(all(centroided(filterMsLevel(x2, 2))))
})


test_that("isCentroidedFromFile", {
    cnt12 <- isCentroidedFromFile(tmt_erwinia_on_disk)
    expect_identical(names(cnt12), featureNames(tmt_erwinia_on_disk))
    expect_identical(isCentroided(tmt_erwinia_on_disk), cnt12)
    cnt1 <- isCentroidedFromFile(tmt_erwinia_on_disk_ms1)
    expect_identical(names(cnt1), featureNames(tmt_erwinia_on_disk_ms1))
    cnt2 <- isCentroidedFromFile(tmt_erwinia_on_disk_ms2)
    expect_identical(names(cnt1), featureNames(tmt_erwinia_on_disk_ms1))
    ##
    ## multiple files
    fls <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_")
    x <- readMSData(fls, mode = "onDisk")
    cnt <- isCentroidedFromFile(x)
    expect_identical(names(cnt), featureNames(x))
    ##
    ## subsetting    
    k <- sort(sample(length(x), 10))
    xx <- x[k]
    expect_identical(isCentroidedFromFile(xx), cnt[k])    
})
