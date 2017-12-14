context("msnset - se conversion methods")

test_that("msnset to se conversion", {
    data(msnset, package = "MSnbase")
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
	data(msnset, package = "MSnbase")
    se <- as(msnset, "SummarizedExperiment")
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
