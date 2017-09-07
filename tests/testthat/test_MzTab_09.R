context("mzTab v0.9")

baseUrl <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/"

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
    expect_warning(prot09 <- readMzTabData(testfile, "PRT", version = "0.9"))
    expect_warning(prot09b <- readMzTabData_v0.9(testfile, "PRT"))
    expect_identical(exprs(prot09), exprs(prot09b))
    expect_identical(fData(prot09), fData(prot09b))    
    expect_warning(pep09 <- readMzTabData(testfile, "PEP", version = "0.9"))
    expect_warning(pep09b <- readMzTabData_v0.9(testfile, "PEP"))
    expect_identical(exprs(pep09), exprs(pep09b))
    expect_identical(fData(pep09), fData(pep09b))
})
