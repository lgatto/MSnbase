context("OnDiskMSnExp class")

inMem <- microtofq_in_mem_ms1
onDisk <- microtofq_on_disk_ms1
multiMsInMem1 <- tmt_im_ms1_sub
multiMsInMem2 <- tmt_im_ms2_sub
multiMsOnDisk <- tmt_od_sub
centroided(inMem) <- TRUE
centroided(onDisk) <- TRUE
centroided(multiMsInMem1) <- TRUE
centroided(multiMsInMem2) <- TRUE
centroided(multiMsOnDisk) <- TRUE

test_that("removePeaks on OnDiskMSnExp with different MS levels", {
    ## o Tests on MSnExp
    multiMsInMem1_rem <- removePeaks(multiMsInMem1)
    expect_true(sum(unlist(intensity(multiMsInMem1_rem)) == 0) >
                sum(unlist(intensity(multiMsInMem1)) == 0))
    multiMsInMem2_rem <- removePeaks(multiMsInMem2)
    expect_true(sum(unlist(intensity(multiMsInMem2_rem)) == 0) >
                sum(unlist(intensity(multiMsInMem2)) == 0))
    ## o Tests on OnDiskMSnExp and comparison with MSnExp.
    multiMsOnDisk_rem <- removePeaks(multiMsOnDisk)
    expect_true(sum(unlist(intensity(multiMsOnDisk_rem)) == 0) >
                sum(unlist(intensity(multiMsOnDisk)) == 0))
    ##   Compare with MSnExp
    expect_true(all.equal(multiMsInMem1_rem,
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 1)))
    expect_true(all.equal(multiMsInMem2_rem,
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 2)))
    ##   Just processing MS 1.
    multiMsOnDisk_rem_1 <- removePeaks(multiMsOnDisk, msLevel. = 1)
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 1),
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 1)))
    spects1 <- spectra(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 2))
    spects2 <- spectra(filterMsLevel(multiMsOnDisk, msLevel. = 2))
    expect_identical(spects1, spects2)
    ## expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_1, msLevel. = 2),
    ##                       filterMsLevel(multiMsOnDisk, msLevel. = 2)))
    ##   Just processing MS 2.
    multiMsOnDisk_rem_2 <- removePeaks(multiMsOnDisk, msLevel. = 2)
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_2, msLevel. = 2),
                          filterMsLevel(multiMsOnDisk_rem, msLevel. = 2)))
    expect_true(all.equal(filterMsLevel(multiMsOnDisk_rem_2, msLevel. = 1),
                          filterMsLevel(multiMsOnDisk, msLevel. = 1)))
})


test_that("clean on OnDiskMSnExp with different MS levels", {
    ## o Tests on MSnExp
    multiMsInMem1_cleaned <- clean(multiMsInMem1)
    expect_true(sum(unlist(intensity(multiMsInMem1_cleaned)) == 0) <
                sum(unlist(intensity(multiMsInMem1)) == 0))
    ## o Tests on OnDiskMSnExp and comparison with MSnExp.
    multiMsOnDisk_cleaned <- clean(multiMsOnDisk)
    expect_true(sum(unlist(intensity(multiMsOnDisk_cleaned)) == 0) <
                sum(unlist(intensity(multiMsOnDisk)) == 0))
    ##   Compare with MSnExp
    expect_true(all.equal(multiMsInMem1_cleaned,
                          filterMsLevel(multiMsOnDisk_cleaned, msLevel. = 1)))
    ##   Just cleaning MS 1.
    multiMsOnDisk_cleaned_1 <- clean(multiMsOnDisk, msLevel. = 1)
    expect_true(all.equal(multiMsOnDisk_cleaned, multiMsOnDisk_cleaned_1))
    ##   Just cleaning MS 2; won't do much at all.
    multiMsOnDisk_cleaned_2 <- clean(multiMsOnDisk, msLevel. = 2)
    expect_true(all.equal(multiMsOnDisk, multiMsOnDisk_cleaned_2))
    ##   Same with msLevel. 4
    multiMsOnDisk_cleaned_4 <- clean(multiMsOnDisk, msLevel. = 4)
    expect_true(all.equal(multiMsOnDisk, multiMsOnDisk_cleaned_4))
})
