test_that(".combineMovingWindow works for Spectrum", {
    ## on a list of spectra.
    spcts <- spectra(tmt_erwinia_in_mem_ms1)
    s_comb <- .combineMovingWindow(spcts)
    expect_equal(length(spcts), length(s_comb))
    expect_equal(unname(lapply(spcts, rtime)), lapply(s_comb, rtime))
    expect_equal(unname(lapply(spcts, msLevel)), lapply(s_comb, msLevel))

    ## Check the first.
    vals_exp <- do.call(rbind, lapply(spcts[1:2], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[1]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[1]]), vals_exp$i)

    ## Check the second.
    vals_exp <- do.call(rbind, lapply(spcts[1:3], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[2]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[2]]), vals_exp$i)

    ## With halfWindowSize 4L
    s_comb <- .combineMovingWindow(spcts, halfWindowSize = 4L)
    expect_equal(length(spcts), length(s_comb))
    expect_equal(unname(lapply(spcts, rtime)), lapply(s_comb, rtime))
    expect_equal(unname(lapply(spcts, msLevel)), lapply(s_comb, msLevel))

    ## Check the first.
    vals_exp <- do.call(rbind, lapply(spcts[1:5], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[1]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[1]]), vals_exp$i)

    ## Check the fifth
    vals_exp <- do.call(rbind, lapply(spcts[1:9], as.data.frame))
    vals_exp <- vals_exp[order(vals_exp$mz), ]
    expect_equal(mz(s_comb[[5]]), vals_exp$mz)
    expect_equal(intensity(s_comb[[5]]), vals_exp$i)
})
