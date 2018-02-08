context("Quantitation with QuantitationParam")


test_that("MS2 isobaric quantitation", {
    data(itraqdata)
    qchar <- quantify(itraqdata, method = "max", reporters = iTRAQ4, qual = FALSE)
    p <- IsobaricTagging(iTRAQ4, method = "max")
    qpar <- quantify(itraqdata, p)    
    qpar@processingData <-
        qchar@processingData <- new("MSnProcess")
    ## all.equal(qchar, qpar)
})

