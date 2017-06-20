test_that("filterEmptySpectra", {
    ## ondisk
    od <- filterFile(microtofq_on_disk_ms1, file = 1)
    od2 <- filterMz(od, mz = c(1005, 1006))
    expect_identical(length(od), length(od2))
    expect_identical(length(filterEmptySpectra(od2)), 0L)
    
    ## inmem
    im <- filterFile(microtofq_in_mem_ms1, file = 1)
    expect_warning(im2 <- filterMz(im, mz = c(1005, 1006)))
    expect_identical(length(im), length(im2))
    expect_identical(length(filterEmptySpectra(im2)), 0L)
})
