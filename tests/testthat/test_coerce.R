data(msnset, package = "MSnbase")

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

context("msnset - se conversion methods")

test_that("msnset to se conversion", {
    se <- as(msnset, "SummarizedExperiment")
    expect_is(se, "SummarizedExperiment")
    expect_equal(colnames(msnset), colnames(se))
    expect_equal(rownames(msnset), rownames(se))
    expect_equal(
        exprs(msnset), 
        SummarizedExperiment::assay(se))
    expect_equal(
        fData(msnset), 
        data.frame(
            SummarizedExperiment::rowData(se), 
            row.names = names(se)))
    expect_equal(
        pData(msnset), 
        data.frame(SummarizedExperiment::colData(se)))
})

test_that("se to msnset conversion", {
    se <- as(msnset, "SummarizedExperiment")
    se <- addMSnSetMetadata(se, msnset)
    msnset_back <- as(se, "MSnSet")
	
    expect_is(msnset_back, "MSnSet")
    expect_equal(colnames(msnset_back), colnames(msnset))
    expect_equal(rownames(msnset_back), rownames(msnset))
    expect_equal(exprs(msnset_back), exprs(msnset))
    expect_equal(fData(msnset_back), fData(msnset))
    expect_equal(pData(msnset_back), pData(msnset))
    
    expect_equal(experimentData(msnset_back), experimentData(msnset))
    expect_equal(protocolData(msnset_back), protocolData(msnset))
    expect_equal(processingData(msnset_back), processingData(msnset))
    expect_equal(qual(msnset_back), qual(msnset))
})
