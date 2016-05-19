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


