test_that(".valid.Backend.processingQueue works", {
    expect_true(is.character(.valid.Backend.processingQueue(list(3, 5))))
    lst <- list(ProcessingStep(mean), ProcessingStep("max"))
    expect_true(is.null(.valid.Backend.processingQueue(lst)))
})
