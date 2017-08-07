test_that("mzRident to data.frame", {
    idf <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid"
    f <- msdata::ident(full.names = TRUE, pattern = idf)
    x <- openIDfile(f)
    iddf <- as(x, "data.frame")
    expect_true(inherits(iddf, "data.frame"))
    expect_true(nrow(iddf) >= length(x))
    iddf2 <- reduce(iddf, key = "spectrumID")
    expect_identical(anyDuplicated(iddf2$spectrumID), 0L)
})

test_that("adding compatible ident with mzID and mzR", {
     quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                      full.name = TRUE, pattern = "mzXML$")
     identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                      full.name = TRUE, pattern = "dummyiTRAQ.mzid")     
     msexp0 <- readMSData(quantFile)
     id1 <- mzR::openIDfile(identFile)
     msexp1 <- addIdentificationData(msexp0, id1)
     id2 <- mzID::mzID(identFile)
     msexp2 <- addIdentificationData(msexp0, id2)

     ## same feature variables are equal
     k <- intersect(fvarLabels(msexp1), fvarLabels(msexp2))     
     expect_equal(fData(msexp1)[, k], fData(msexp2)[, k])

     expect_equal(fData(msexp1)[, "chargeState"], fData(msexp2)[, "chargestate"])
     expect_equal(fData(msexp1)[, "passThreshold"], fData(msexp2)[, "passthreshold"])
     expect_equal(fData(msexp1)[, "experimentalMassToCharge"],
                  fData(msexp2)[, "experimentalmasstocharge"])
     expect_equal(fData(msexp1)[, "calculatedMassToCharge"],
                  fData(msexp2)[, "calculatedmasstocharge"])
     expect_equal(fData(msexp1)[, "sequence"], fData(msexp2)[, "pepseq"])
     expect_equal(fData(msexp1)[, "DatabaseAccess"],
                  fData(msexp2)[, "accession"])
     ## expect_equal(fData(msexp1)[, "DatabaseDescription"],
     ##              fData(msexp2)[, "description"]) ## checked manually
     expect_equal(fData(msexp1)[, "MS.GF.EValue"],
                  fData(msexp2)[, "ms.gf.evalue"])
     expect_equal(fData(msexp1)[, "MS.GF.DeNovoScore"],
                  fData(msexp2)[, "ms.gf.denovoscore"])
     expect_equal(fData(msexp1)[, "MS.GF.RawScore"],
                  fData(msexp2)[, "ms.gf.rawscore"])
     expect_equal(fData(msexp1)[, "MS.GF.SpecEValue"],
                  fData(msexp2)[, "ms.gf.specevalue"])
     expect_equal(fData(msexp1)[, "isDecoy"], fData(msexp2)[, "isdecoy"])
     expect_equal(fData(msexp1)[, "DBseqLength"], fData(msexp2)[, "length"])     
})

test_that("adding compatible ident to MSnExp and MSnSet", {
     quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                      full.name = TRUE, pattern = "mzXML$")
     identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                      full.name = TRUE, pattern = "dummyiTRAQ.mzid")
     msexp <- readMSData2(quantFile)
     msset <- quantify(msexp, method = "max", reporters = iTRAQ4)
     msexp <- addIdentificationData(msexp, identFile)
     msset <- addIdentificationData(msset, identFile) 
})
