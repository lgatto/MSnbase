context("iPQF feature aggregation")

test_that("data is valid", {
    ## test input
    data("msnset2")
    expect_true(validObject(msnset2))
    ## test output
    fl <- list.files(system.file("extdata", package = "MSnbase"),
                     full.names=TRUE, pattern = "ipqfres.rda")
    nms <- load(fl)
    val <- sapply(nms, function(x) validObject(get(x)))
    expect_true(all(val))
})

test_that("iPQF unit test", {
    data("msnset2")
    fl <- list.files(system.file("extdata", package = "MSnbase"),
                     full.names = TRUE, pattern = "ipqfres.rda")
    load(fl) ## ipqf1 to 5

    res1 <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            fun = "iPQF",
                            low.support.filter = FALSE,
                            ratio.calc = "sum",
                            cv = TRUE,
                            method.combine = FALSE)
    res1@processingData <-
        ipqf1@processingData <- new("MSnProcess")
    res1@experimentData <-
        ipqf1@experimentData <- new("MIAPE")
    res1@.__classVersion__ <- ipqf1@.__classVersion__
    res1 <- droplevels(res1)   ## new
    ipqf1 <- droplevels(ipqf1) ## new
    expect_equal(res1, ipqf1)

    dflt <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            cv = TRUE,
                            fun = "iPQF")
    dflt@processingData <-
        res1@processingData <- new("MSnProcess")
    dflt@experimentData <-
        res1@experimentData <- new("MIAPE")
    dflt@.__classVersion__ <- res1@.__classVersion__
    dflt <- droplevels(dflt)
    res1 <- droplevels(res1)
    expect_equal(res1, dflt)
    
    res2 <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            fun = "iPQF",
                            cv = TRUE,
                            low.support.filter = FALSE,
                            ratio.calc = "sum",
                            method.combine = TRUE)
    res2@processingData <-
        ipqf2@processingData <- new("MSnProcess")
    res2@experimentData <-
        ipqf2@experimentData <- new("MIAPE")
    res2@.__classVersion__ <- ipqf2@.__classVersion__
    res2 <- droplevels(res2)   ## new
    ipqf2 <- droplevels(ipqf2) ## new
    expect_equal(res2, ipqf2)

    res3 <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            fun = "iPQF",
                            cv = TRUE,
                            low.support.filter = TRUE,
                            ratio.calc = "sum",
                            method.combine = FALSE)
    res3@processingData <-
        ipqf3@processingData <- new("MSnProcess")
    res3@experimentData <-
        ipqf3@experimentData <- new("MIAPE")
    res3@.__classVersion__ <- ipqf3@.__classVersion__
    res3 <- droplevels(res3)
    ipqf3 <- droplevels(ipqf3)
    expect_equal(res3, ipqf3)

    res4 <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            fun = "iPQF",
                            cv = TRUE,
                            low.support.filter = TRUE,
                            ratio.calc = "none",
                            method.combine = FALSE)
    res4@processingData <-
        ipqf4@processingData <- new("MSnProcess")
    res4@experimentData <-
        experimentData(ipqf4) <- new("MIAPE")
    res4@.__classVersion__ <- ipqf4@.__classVersion__
    res4 <- droplevels(res4)
    ipqf4 <- droplevels(ipqf4)
    expect_equal(res4, ipqf4)

    res5 <- combineFeatures(msnset2,
                            groupBy = fData(msnset2)$accession,
                            redundancy.handler = "unique",
                            fun = "iPQF",
                            cv = TRUE,
                            low.support.filter = TRUE,
                            ratio.calc = "X114.ions",
                            method.combine = FALSE)
    res5@processingData <-
        ipqf5@processingData <- new("MSnProcess")
    res5@experimentData <-
        ipqf5@experimentData <- new("MIAPE")
    res5@.__classVersion__ <- ipqf5@.__classVersion__
    res5 <- droplevels(res5)
    ipqf5 <- droplevels(ipqf5)
    expect_equal(res5, ipqf5)    
})
