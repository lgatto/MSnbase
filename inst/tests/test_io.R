context("IO testing")

test_that("readMSData >1 files", {
  file1 <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  tdir <- tempdir()
  file2 <- file.path(tdir, sub("iTRAQ", "iTRAQ2", basename(file1)))
  stopifnot(file.copy(file1, file2, overwrite = TRUE))
  raw <- readMSData(c(file1, file2))
  raw
})


test_that("Compatibility between writeMgfData and readMgfData", {
  data(itraqdata)
  tf <- tempfile()
  ## no spectra order issue here, as < 10 spectra
  ## see [read|write]MgfData for details
  d1 <- itraqdata[1:3]
  writeMgfData(d1, con = tf)
  d2 <- readMgfData(tf, verbose = FALSE)
  expect_true(all(sapply(1:3,
                         function(i) all.equal(as.data.frame(d1[[i]]),
                                               as.data.frame(d1[[i]])))))
  expect_true(all.equal(precursorMz(d1), precursorMz(d1)))
  unlink(tf)
})

## Removing this test, as readMzXML has been defunct
## in MSnbase version 1.3.5
## test_that("Compare readMzXMLData and readMSData output", {
##   file <- dir(system.file(package="MSnbase", dir="extdata"),
##               full.names=TRUE,
##               pattern="mzXML$")
##   expect_warning(aa <- readMzXMLData(file, verbose=FALSE))
##   bb <- readMSData(file, verbose=FALSE)  
##   ## comparing data
##   expect_true(all.equal(header(aa), header(bb)))
##   expect_true(all.equal(intensity(aa), intensity(bb)))
##   expect_true(all.equal(mz(aa), mz(bb)))  
## })

test_that("Testing write.exprs and readMSnSet", {
  data(itraqdata)
  colchars <- c("ProteinAccession", "PeptideSequence", "retention.time", "precursor.mz")
  tf <- tempfile()
  x <- quantify(itraqdata, reporters = iTRAQ4, method = "max", verbose = FALSE)
  write.exprs(x, file = tf)  
  y <- readMSnSet(tf)
  expect_true(all.equal(exprs(x), exprs(y)))
  ## unlink(tf)
  write.exprs(x, fDataCols = colchars, file = tf)
  tmp <- read.table(tf)
  expect_true(all(dim(tmp) == c(nrow(x), ncol(x) + length(colchars))))
  expect_true(all(colnames(tmp) == c(sampleNames(x), colchars)))
  expect_true(all(rownames(tmp) == featureNames(x)))
  unlink(tf)
})
