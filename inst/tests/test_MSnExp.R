context("MSnExp class")

test_that("readMzXMLData and dummy MSnExp instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  expect_that(class(aa)=="MSnExp",is_true())
  expect_that(length(aa),equals(98))
  ## MS levels
  expect_that(length(msLevel(aa)),equals(98))
  expect_that(unique(msLevel(aa)),equals(2))
  expect_that(length(ms1scan(aa)),equals(98))
  expect_that(length(unique(ms1scan(aa))),equals(7))
  ## Precursor MZ
  expect_that(length(precursorMz(aa)),equals(98))
  expect_that(precursorMz(aa)[1],is_a("numeric"))
  expect_that(length(unique(precursorMz(aa))),equals(15))
  expect_that(range(precursorMz(aa)),
              equals(c(424.76651001,1007.04766850)))
  expect_that(as.numeric(sort(precursorMz(aa))[1]),equals(424.76651001)) ## [*]
  ## Retention time
  expect_that(length(rtime(aa)),equals(98))
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_that(range(rtime(aa)),
              equals(c(1982.92,6179.02)))
  expect_that(as.numeric(sort(rtime(aa))[1]),equals(1982.92)) ## [*]
  ## [*] using as.numeric because rtime and precursorMz return named numerics
})


context("data integrity")

test_that("spectra order and integrity", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  sp <- assayData(aa)$X98
  nms.0 <- ls(assayData(aa))
  nms.1 <- ls(assayData(removePeaks(aa,verbose=FALSE)))
  expect_that(nms.1,equals(nms.0))
  expect_that(all.equal(removePeaks(sp,t=2),
            assayData(removePeaks(aa,t=2,verbose=FALSE))$X98),
            is_true())
  expect_that(tic(sp),equals(334))
  expect_that(tic(removePeaks(sp,1)),equals(56))
  expect_that(tic(removePeaks(sp,2)),equals(22))
  expect_that(tic(removePeaks(sp,3)),equals(0))
  expect_that(peaksCount(sp),equals(849))
  expect_that(peaksCount(clean(removePeaks(sp,1))),equals(65))
  expect_that(peaksCount(clean(removePeaks(sp,2))),equals(21))
  expect_that(peaksCount(clean(removePeaks(sp,3))),equals(0))
  expect_that(tic(clean(removePeaks(sp,0))),equals(tic(sp)))
})

context("quantification and MSnSet instance")


test_that("quantification", {
  ## dummy Spectrum
  int <- c(0,2,3,1,0)
  mz <- c(114.05,
          114.10,
          114.13,
          114.15,
          114.16)
  sp <- new("Spectrum2",
            intensity=int,
            mz=mz)
  expect_that(validObject(sp),is_true())
  data(iTRAQ4)
  ##quantify(sp,iTRAQ4[1],"sum")
})
