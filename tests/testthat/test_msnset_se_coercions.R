context("msnset - se conversion methods")

test_that("msnset to se conversion", {
	data(msnset, package = "MSnbase")
	se <- as(msnset, "SummarizedExperiment")
	
	expect_is(se, "SummarizedExperiment")
	expect_equal(colnames(msnset), colnames(se))
	expect_equal(rownames(msnset), rownames(se))
	expect_equal(exprs(msnset), assay(se))
	expect_equal(fData(msnset), data.frame(rowData(se), row.names = names(se)))
	expect_equal(pData(msnset), data.frame(colData(se)))
})

test_that("se to msnset conversion", {
	data(msnset, package = "MSnbase")
    se <- as(msnset, "SummarizedExperiment")
	msnset_back <- as(se, "MSnSet")
	
	expect_is(msnset_back, "MSnSet")
	expect_equal(colnames(msnset_back), colnames(se))
	expect_equal(rownames(msnset_back), rownames(se))
	expect_equal(exprs(msnset_back), assay(se))
	expect_equal(fData(msnset_back), data.frame(rowData(se), row.names = names(se)))
	expect_equal(pData(msnset_back), data.frame(colData(se)))
})
