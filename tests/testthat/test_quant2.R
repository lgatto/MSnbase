context("Quantitation with QuantitationParam")


test_that("MS2 isobaric quantitation MSnExp", {
    data(itraqdata)
    qchar <- quantify(itraqdata, method = "max", reporters = iTRAQ4, qual = FALSE)
    p <- IsobaricTagging(iTRAQ4, method = "max")
    qpar <- quantify(itraqdata, p)    
    qpar@processingData <-
        qchar@processingData <- new("MSnProcess")
    expect_equivalent(qchar, qpar)
})

test_that("MS2 isobaric quantitation OnDiskMSnExp", {
     f <- dir(system.file(package = "MSnbase", dir = "extdata"),
              full.name = TRUE,
              pattern = "mzXML$")
     x <- readMSData(f, mode = "onDisk")
     p <- IsobaricTagging(iTRAQ4, method = "max")     
     qpar <- quantify(x, p)
     qchar <- quantify(x, method = "max", reporters = iTRAQ4)
     qchar@processingData <-
         qpar@processingData <- new("MSnProcess")
     ## remove on feature variable
     fData(qchar)$reporterMzs <- NULL
     expect_equal(qpar, qchar)
})


test_that("MSnExp vs OnDiskMSnExp", {
     f <- dir(system.file(package = "MSnbase", dir = "extdata"),
              full.name = TRUE,
              pattern = "mzXML$")
     od <- readMSData(f, mode = "onDisk", msLevel = 2L)
     im <- readMSData(f, mode = "inMemory")
     p <- IsobaricTagging(iTRAQ4, method = "max")
     qod <- quantify(od, p)
     qim <- quantify(im, p) 
     ## these two have different processingData (expected) but also
     ## different fData - od has many more feature variables
     expect_identical(exprs(qim), exprs(qod))
})


