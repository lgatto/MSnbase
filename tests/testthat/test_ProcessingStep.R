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

