context("msnset - se conversion methods")

test_that("msnset2se conversion", {
	data(msnset, package = "MSnbase")
	se <- msnset2se(msnset)
	
	expect_is(se, "SummarizedExperiment")
	expect_equal(colnames(msnset), colnames(se))
	expect_equal(rownames(msnset), rownames(se))
	expect_equal(Biobase::exprs(msnset), SummarizedExperiment::assay(se))
	expect_equal(Biobase::fData(msnset), data.frame(SummarizedExperiment::rowData(se), row.names = se@NAMES))
	expect_equal(Biobase::pData(msnset), data.frame(SummarizedExperiment::colData(se)))
})

test_that("se2msnset conversion", {
	data(msnset, package = "MSnbase")
	se <- msnset2se(msnset)
	msnset_back <- se2msnset(se)
	
	expect_is(msnset_back, "MSnSet")
	expect_equal(colnames(msnset_back), colnames(se))
	expect_equal(rownames(msnset_back), rownames(se))
	expect_equal(Biobase::exprs(msnset_back), SummarizedExperiment::assay(se))
	expect_equal(Biobase::fData(msnset_back), data.frame(SummarizedExperiment::rowData(se), row.names = se@NAMES))
	expect_equal(Biobase::pData(msnset_back), data.frame(SummarizedExperiment::colData(se)))
})