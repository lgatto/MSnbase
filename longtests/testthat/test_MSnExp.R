test_that("estimateMzScattering works", {
    expect_error(estimateMzScattering(4))

    res <- estimateMzScattering(tmt_erwinia_in_mem_ms1)
    mzr <- estimateMzResolution(tmt_erwinia_in_mem_ms1)
    idx <- which.max(spectrapply(tmt_erwinia_in_mem_ms1, peaksCount))
    ## m/z scattering should be smaller than m/z resolution
    expect_true(res[[idx]] < mzr[[idx]])

    ## .estimate_mz_scattering_list.
    res_2 <- .estimate_mz_scattering_list(spectra(tmt_erwinia_in_mem_ms1),
                                          timeDomain = FALSE)
    expect_equal(unname(res), res_2)
})
