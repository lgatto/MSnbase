test_that("All defunct functions", {
    expect_error(readMzXMLData(), "Defunct")
    expect_error(readMzXMLData(), "Defunct")
    expect_error(makeMTD(), "Defunct")
    expect_error(makePRT(), "Defunct")
    expect_error(makePEP(), "Defunct")
    expect_error(writeMzTabData(), "Writing support for mzTab is defunct.")
    expect_error(writeMzTabData(), "Writing support for mzTab is defunct.")
})
