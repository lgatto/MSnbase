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
     expect_equal(qpar, qchar)
})
