data(msnset)

test_that("coerce experimentData", {
    x <- as(experimentData(msnset), "list")
    expect_length(x, length(slotNames(experimentData(msnset))) - 1)
    x <- MSnbase:::.reduce_list(x)
    expect_length(x, 5)
})


test_that("coerce protocolData", {
    x <- as(protocolData(msnset), "list")
    expect_length(x, 0)
})


test_that("coerce MSnProcess", {
    x <- as(processingData(msnset), "list")
    expect_length(x, length(slotNames(processingData(msnset))) - 1)
    x <- MSnbase:::.reduce_list(x)
    expect_length(x, 4)
})
