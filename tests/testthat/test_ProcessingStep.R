############################################################
## Testing ProcessingStep functionality.
test_that("ProcessingStep constructor", {
    ps <- new("ProcessingStep", FUN="mean", ARGS=list(na.rm=TRUE))
    expect_that(new("ProcessingStep", FUN="aaaa"),
                throws_error("Function 'aaaa' not found!"))
    ## aaaa <- 1
    ## expect_that(new("ProcessingStep", FUN="aaaa"),
    ##             throws_error("'aaaa' is not a function!"))

    ProcessingStep("mean", list(na.rm=TRUE))
})

test_that("ProcessingStep executeProcessingStep", {
    ps <- ProcessingStep("sum", list(c(1, 4, 5)))
    expect_identical(executeProcessingStep(ps), 10)

    ps <- ProcessingStep("sum", list(c(1, 4, 5, NA)))
    ## Pass optional arguments to ...
    expect_identical(executeProcessingStep(ps, na.rm = TRUE), 10)
})

test_that(".apply_processing_queue works", {
    sps <- spectra(sciex)
    expect_equal(sps, MSnbase:::.apply_processing_queue(sps))
    q <- list(ProcessingStep("fromFile"))
    res <- MSnbase:::.apply_processing_queue(sps, q)
    expect_equal(res, lapply(sps, fromFile))

    q <- list(ProcessingStep(FUN = function(x, a) {
        fromFile(x) * a
    }, ARGS = list(a = 4)))
    res <- MSnbase:::.apply_processing_queue(sps, q)
    expect_equal(unlist(res), unlist(lapply(sps, fromFile)) * 4)

    q <- list(ProcessingStep("removePeaks", list(t = 0)),
              ProcessingStep("clean", list(all = TRUE)),
              ProcessingStep("intensity"))
    res <- MSnbase:::.apply_processing_queue(sps[1:3], q)
    expect_equal(res[[1]], intensity(clean(removePeaks(sps[[1]], t = 20),
                                           all = TRUE)))
})
