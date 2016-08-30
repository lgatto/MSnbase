context("readWriteMgfData")

test_that("extractMgfSpectrum2Info", {
  mgf <- c("TITLE=foobar File=\"foobar.raw\", Native ID:\"controllerType=0 controllerNumber=1, scan=100\"",
           "RTINSECONDS=600",
           "PEPMASS=100 50000",
           "CHARGE=3+",
           "10 100",
           "11 200",
           "12 300",
           "13 400")

  s <- new("Spectrum2",
           rt = 600,
           scanIndex = 0,
           precursorMz = 100,
           precursorIntensity = 50000,
           precursorCharge = 3,
           mz = 10:13,
           intensity = seq(100, 400, by=100),
           fromFile = 1L,
           centroided = TRUE)
  fdata <- c(TITLE="foobar File=\"foobar.raw\", Native ID:\"controllerType=0 controllerNumber=1, scan=100\"",
             RTINSECONDS="600",
             PEPMASS="100 50000",
             CHARGE="3+")

  result <- list(spectrum=s, fdata=fdata)

  expect_identical(MSnbase:::extractMgfSpectrum2Info(mgf, centroided = TRUE), result)
})
