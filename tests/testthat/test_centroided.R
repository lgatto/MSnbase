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
    ## multiple files; use the microtofq files
    cnt <- isCentroidedFromFile(microtofq_on_disk)
    expect_identical(names(cnt), featureNames(microtofq_on_disk))
    ##
    ## subsetting
    set.seed(123) ## see issue #338
    k <- sort(sample(length(microtofq_on_disk), 10))
    xx <- microtofq_on_disk[k]
    expect_identical(isCentroidedFromFile(xx), cnt[k])
})
