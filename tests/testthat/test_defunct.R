test_that("All defunct functions", {
    expect_error(.Defunct(readMzXMLData()), "Defunct")
    expect_error(.Defunct(readMzXMLData()), "Defunct")
    expect_error(.Defunct(makeMTD()), "Defunct")
    expect_error(.Defunct(makePRT()), "Defunct")
    expect_error(.Defunct(makePEP()), "Defunct")
    expect_error(.Defunct(writeMzTabData()), "Defunct")
})
