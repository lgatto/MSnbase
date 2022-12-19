context("MzTab class")

baseUrl <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/1_0-Proteomics-Release/"
baseUrlM <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/2_0-Metabolomics_Release/"

test_that("MzTab creation and accessors", {
    fl <- "iTRAQ_CQI.mzTab"
    fl <- file.path(baseUrl, fl)
    xx <- MzTab(fl)
    expect_true(validObject(xx))
    expect_null(show(xx))
    ## Accessors
    expect_is(metadata(xx), "list")
    expect_identical(length(metadata(xx)), 65L)
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

test_that("MzTab-M creation and accessors", {
    fl <- "MTBLS263.mztab"
    fl <- file.path(baseUrlM, fl)
    xx <- MzTab(fl)
    expect_true(validObject(xx))
    expect_null(show(xx))
    ## Accessors
    expect_is(metadata(xx), "list")
    expect_identical(length(metadata(xx)), 74L)

    ## In mzTab-M 2.0 metadata(xx)$`mzTab-mode`=="NULL"
    expect_identical(mzTabMode(xx), metadata(xx)$`mzTab-mode`)

    ## In mzTab-M 2.0 metadata(xx)$`mzTab-type`=="NULL"
    expect_identical(mzTabType(xx), metadata(xx)$`mzTab-type`)

    expect_identical(fileName(xx), fl)
    expect_is(proteins(xx), "data.frame")
    expect_identical(dim(proteins(xx)), c(0L, 0L))

    expect_is(peptides(xx), "data.frame")
    expect_identical(dim(peptides(xx)), c(0L, 0L))

    expect_is(psms(xx), "data.frame")
    expect_identical(dim(psms(xx)), c(0L, 0L))

    expect_is(smallMolecules(xx), "data.frame")
    expect_identical(dim(smallMolecules(xx)), c(17L, 24L))
    expect_is(moleculeFeatures(xx), "data.frame")
    expect_identical(dim(moleculeFeatures(xx)), c(19L, 24L))
    expect_is(moleculeEvidence(xx), "data.frame")
    expect_identical(dim(moleculeEvidence(xx)), c(19L, 24L))

    expect_is(comments(xx), "character")
    expect_identical(length(comments(xx)), 0L)
})

test_that("MzTab-M MS-Dial import", {
  fl <- "msdial/lcmsms_swath_lipid_height_mzTab.mztab"
  fl <- file.path(baseUrlM, fl)
  xx <- MzTab(fl)
  expect_true(validObject(xx))
  expect_null(show(xx))

  ## Accessors
  expect_is(metadata(xx), "list")
  expect_identical(length(metadata(xx)), 178L)

  ## In mzTab-M 2.0 metadata(xx)$`mzTab-mode`=="NULL"
  expect_identical(mzTabMode(xx), metadata(xx)$`mzTab-mode`)

  ## In mzTab-M 2.0 metadata(xx)$`mzTab-type`=="NULL"
  expect_identical(mzTabType(xx), metadata(xx)$`mzTab-type`)

  expect_identical(fileName(xx), fl)
  expect_is(proteins(xx), "data.frame")
  expect_identical(dim(proteins(xx)), c(0L, 0L))

  expect_is(peptides(xx), "data.frame")
  expect_identical(dim(peptides(xx)), c(0L, 0L))

  expect_is(psms(xx), "data.frame")
  expect_identical(dim(psms(xx)), c(0L, 0L))

  expect_is(smallMolecules(xx), "data.frame")
  expect_identical(dim(smallMolecules(xx)), c(1172L, 47L))
  expect_is(moleculeFeatures(xx), "data.frame")
  expect_identical(dim(moleculeFeatures(xx)), c(1172L, 34L))
  expect_is(moleculeEvidence(xx), "data.frame")
  expect_identical(dim(moleculeEvidence(xx)), c(199L, 23L))

  expect_is(comments(xx), "character")
  expect_identical(length(comments(xx)), 1L)
})


test_that("MzTab reading and writing", {
    fl <- "iTRAQ_CQI.mzTab"
    fl <- file.path(baseUrl, fl)
    xx <- MzTab(fl)

    outfile <- tempfile()
    writeMzTabData(xx, file = outfile)

    inlines <- readLines(fl)
    inlines <- inlines[!grepl("^COM", inlines)]
    outlines <- readLines(outfile)

    #    expect_true(all.equal(inlines, outlines))
})

test_that("MzTab-M reading and writing", {
    fl <- "MTBLS263.mztab"
    fl <- file.path(baseUrlM, fl)
    xx <- MzTab(fl)

    outfile <- tempfile()
    writeMzTabData(xx, file = outfile)

    inlines <- readLines(fl)
    inlines <- inlines[!grepl("^COM", inlines)]
    outlines <- readLines(outfile)

    #    expect_true(all.equal(inlines, outlines))
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
