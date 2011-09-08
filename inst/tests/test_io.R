context("IO testing")

test_that("Compatibility between writeMgfData and readMgfData", {
  data(itraqdata)
  tf <- tempfile()
  ## no spectra order issue here, as < 10 spectra
  ## see [read|write]MgfData for details
  d1 <- itraqdata[1:3]
  writeMgfData(d1,file=tf)
  d2 <- readMgfData(tf, verbose=FALSE)
  expect_true(all(sapply(1:3, function(i) all.equal(as.data.frame(d1[[i]]),as.data.frame(d1[[i]])))))
  expect_true(all.equal(precursorMz(d1),precursorMz(d1)))
  unlink(tf)
})

test_that("Compare readMzXMLData and readMSData output", {
  file <- dir(system.file(package="MSnbase", dir="extdata"),
              full.name=TRUE,
              pattern="mzXML$")
  aa <- readMzXMLData(file, verbose=FALSE)  
  bb <- readMSData(file, verbose=FALSE)  
  ## comparing data
  expect_true(all.equal(header(aa),header(bb)))
  expect_true(all.equal(intensity(aa),intensity(bb)))
  expect_true(all.equal(mz(aa),mz(bb)))  
})
