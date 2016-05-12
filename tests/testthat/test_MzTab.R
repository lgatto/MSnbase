context("MzTab class")

baseUrl <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/"

test_that("MzTab creation and accessors", {
    fl <- "iTRAQ_CQI.mzTab"
    fl <- file.path(baseUrl, fl)
    xx <- MzTab(fl)
    expect_true(validObject(xx))
    expect_null(show(xx))
    ## Accessors
    expect_is(metadata(xx), "list")
    expect_identical(length(metadata(xx)), 64L)
    expect_identical(mzTabMode(xx), metadata(xx)$`mzTab-mode`)
    expect_identical(mzTabMode(xx), "Complete")
    expect_identical(mzTabType(xx), metadata(xx)$`mzTab-type`)
    expect_identical(mzTabType(xx), "Quantification")
    expect_identical(fileName(xx), fl)
    expect_is(proteins(xx), "data.frame")
    expect_identical(dim(proteins(xx)), c(5L, 55L))
    expect_is(peptides(xx), "data.frame")
    expect_identical(dim(peptides(xx)), c(0L, 0L))
    expect_is(psms(xx), "data.frame")
    expect_identical(dim(psms(xx)), c(36L, 18L))
    expect_is(smallMolecules(xx), "data.frame")
    expect_identical(dim(smallMolecules(xx)), c(0L, 0L))
    expect_is(comments(xx), "character")
    expect_identical(length(comments(xx)), 5L)
})

test_that("Conversion to MSnSetList", {
    fl <- "iTRAQ_CQI.mzTab"
    fl <- file.path(baseUrl, fl)
    xx <- MzTab(fl)
    msl <- as(xx, "MSnSetList")
    expect_true(validObject(msl))
    expect_null(show(msl))
    expect_true(length(msl) == 3L)
    dims0 <- structure(list(Proteins = c(5L, 16L),
                            Peptides = c(0L, 0L),
                            PSMs = c(36L, 0L)))
    expect_identical(lapply(msl, dim), dims0)
    ##
    ms <- msl[['Proteins']]
    dfr <- proteins(xx)
    expect_identical(ncol(dfr), length(fvarLabels(ms)) + ncol(ms))
    ms <- msl[['Peptides']]
    dfr <- peptides(xx)
    expect_identical(ncol(dfr), length(fvarLabels(ms)) + ncol(ms))
    ms <- msl[['PSMs']]
    dfr <- psms(xx)
    expect_identical(ncol(dfr), length(fvarLabels(ms)) + ncol(ms))              
})

test_that("read MzTab data v 0.9 and 1.0", {
    f <- "PRIDE_Exp_Complete_Ac_16649.xml-mztab.txt"
    testfile <- file.path(baseUrl, f)
    prot <- readMzTabData(testfile, "PRT")
    psm <- readMzTabData(testfile, "PSM")
    mzt <- MzTab(testfile)
    ml <- as(mzt, "MSnSetList")
    expect_identical(exprs(ml[["Proteins"]]), exprs(prot))
    expect_identical(exprs(ml[["PSMs"]]), exprs(psm))
    expect_identical(fData(ml[["Proteins"]]), fData(prot))
    expect_identical(fData(ml[["PSMs"]]), fData(psm))

    testfile <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/legacy/jmztab-1.0/examples/mztab_itraq_example.txt"
    expect_error(readMzTabData(testfile, "PSM", version = "0.9"))
    prot09 <- readMzTabData(testfile, "PRT", version = "0.9")
    prot09b <- readMzTabData_v0.9(testfile, "PRT")
    expect_identical(exprs(prot09), exprs(prot09b))
    expect_identical(fData(prot09), fData(prot09b))    
    pep09 <- readMzTabData(testfile, "PEP", version = "0.9")
    pep09b <- readMzTabData_v0.9(testfile, "PEP")
    expect_identical(exprs(pep09), exprs(pep09b))
    expect_identical(fData(pep09), fData(pep09b))
})
