context("[OnDisk]MSnExp equality")

test_that("Equality function", {
    ## testing it on spectra only.
    expect_true(all.equal(tmt_erwinia_in_mem_ms1, tmt_erwinia_on_disk_ms1))
    expect_true(all.equal(tmt_erwinia_in_mem_ms2, tmt_erwinia_on_disk_ms2))
    ## postive controls
    expect_true(all.equal(tmt_erwinia_in_mem_ms1, tmt_erwinia_in_mem_ms1))
    expect_true(all.equal(tmt_erwinia_in_mem_ms2, tmt_erwinia_in_mem_ms2))
    expect_true(all.equal(tmt_erwinia_on_disk_ms1, tmt_erwinia_on_disk_ms1))
    ## negative controls
    expect_false(isTRUE(all.equal(tmt_erwinia_in_mem_ms1,
                                  tmt_erwinia_in_mem_ms2)))
    expect_false(isTRUE(all.equal(tmt_erwinia_on_disk_ms1,
                                  tmt_erwinia_on_disk_ms2)))
})
