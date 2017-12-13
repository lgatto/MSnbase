context("fdata-selection")

test_that("Feature variable selection", {
    data(hyperLOPIT2015, package = "pRolocdata")
    fv <- fvarLabels(hyperLOPIT2015)
    i <- sort(sample(length(fv), 10))
    k <- fv[i]
    l <- logical(length(fv))
    l[i] <- TRUE
    expect_equal(selectFeatureData(hyperLOPIT2015, fcol = i),
                 selectFeatureData(hyperLOPIT2015, fcol = k))
    expect_equal(selectFeatureData(hyperLOPIT2015, fcol = i),
                 selectFeatureData(hyperLOPIT2015, fcol = l))
})

test_that("selectFeatureData,OnDiskMSnExp works", {
    on_disk <- tmt_erwinia_on_disk
    fvars <- c("fileIdx", "spIdx", "acquisitionNum", "retentionTime", "msLevel",
               "precursorScanNum", "centroided")
    res <- selectFeatureData(on_disk, fcol = fvars)
    expect_equal(fvarLabels(res), fvars)
    expect_true(validObject(res))
    sp_full <- on_disk[[3]]
    sp_red <- res[[3]]
    expect_equal(mz(sp_full), mz(sp_red))
    expect_equal(intensity(sp_full), intensity(sp_red))
    expect_equal(rtime(sp_full), rtime(sp_red))
    expect_equal(acquisitionNum(sp_full), acquisitionNum(sp_red))
    expect_equal(msLevel(sp_full), msLevel(sp_red))
    expect_equal(centroided(sp_full), centroided(sp_red))
    expect_equal(precScanNum(sp_full), precScanNum(sp_red))
    ## Stuff that is no longer set.
    expect_true(is.na(tic(sp_red)))
    expect_true(!is.na(tic(sp_full)))
    expect_true(is.na(polarity(sp_red)))
    expect_true(!is.na(polarity(sp_full)))

    ## And for MS1 data
    on_disk <- microtofq_on_disk
    res <- selectFeatureData(on_disk, fcol = fvars)
    expect_equal(fvarLabels(res), fvars)
    expect_true(validObject(res))
    sp_full <- on_disk[[3]]
    sp_red <- res[[3]]
    expect_equal(mz(sp_full), mz(sp_red))
    expect_equal(intensity(sp_full), intensity(sp_red))
    expect_equal(rtime(sp_full), rtime(sp_red))
    expect_equal(acquisitionNum(sp_full), acquisitionNum(sp_red))
    expect_equal(msLevel(sp_full), msLevel(sp_red))
    expect_equal(centroided(sp_full), centroided(sp_red))
    ## Stuff that is no longer set.
    expect_true(is.na(tic(sp_red)))
    expect_true(!is.na(tic(sp_full)))
    expect_true(is.na(polarity(sp_red)))
    expect_true(!is.na(polarity(sp_full)))
    sp_full <- on_disk[[2]]
    sp_red <- res[[2]]
    expect_equal(mz(sp_full), mz(sp_red))
    expect_equal(intensity(sp_full), intensity(sp_red))
    expect_equal(rtime(sp_full), rtime(sp_red))
    expect_equal(acquisitionNum(sp_full), acquisitionNum(sp_red))
    expect_equal(msLevel(sp_full), msLevel(sp_red))
    expect_equal(centroided(sp_full), centroided(sp_red))
    ## Stuff that is no longer set.
    expect_true(is.na(tic(sp_red)))
    expect_true(!is.na(tic(sp_full)))
    expect_true(is.na(polarity(sp_red)))
    expect_true(!is.na(polarity(sp_full)))
})

test_that("selectFeatureData,MSnExp works", {
    in_mem <- tmt_erwinia_in_mem_ms1
    ## fData is currently a single column data.frame...
})

test_that("requiredFvarLabels works", {
    expect_error(requiredFvarLabels(x = "other"))
    res <- requiredFvarLabels("OnDiskMSnExp")
    expect_equal(res, MSnbase:::.MSnExpReqFvarLabels)
    expect_equal(requiredFvarLabels("MSnExp"), character())
    expect_equal(requiredFvarLabels("MSnSet"), character())
})
