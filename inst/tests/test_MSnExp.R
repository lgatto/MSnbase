context("MSnExp class")

test_that("MSnExp validity", {
  expect_true(validObject(new("MSnExp")))
  })

test_that("readMzXMLData and dummy MSnExp msLevel 2 instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  ## aa <- readMzXMLData(file,verbose=FALSE)
  aa <- readMSData(file,verbose=FALSE)
  expect_true(class(aa)=="MSnExp")
  ## centroided get and set
  expect_false(any(centroided(aa)))
  val <- rep(TRUE,length(aa))
  centroided(aa) <- val
  expect_true(validObject(aa))
  expect_true(all(centroided(aa)))
  val[sample(length(aa),2)] <- FALSE
  centroided(aa) <- val
  expect_true(sum(centroided(aa))==length(aa)-2)
  centroided(aa) <- rep(FALSE,length(aa))
  expect_false(any(centroided(aa)))  
  ## checking slots and methods
  expect_equal(length(aa),5)
  expect_that(nrow(header(aa)),equals(length(aa)))
  expect_that(names(header(aa)),
              equals(c("index","file","retention.time",
                       "precursor.mz","peaks.count","tic",
                       "ms.level","charge","collision.energy")))
  ## MS levels
  expect_equal(length(msLevel(aa)),5)
  expect_equal(unique(msLevel(aa)),2)
  expect_equal(length(precScanNum(aa)),5)
  expect_equal(length(unique(precScanNum(aa))),1)
  ## Precursor MZ
  expect_equal(length(precursorMz(aa)),5)
  expect_that(precursorMz(aa)[1],is_a("numeric"))
  expect_equal(length(unique(precursorMz(aa))),4)
  expect_equal(range(precursorMz(aa)),
               c(437.80401611,716.34051514))
  expect_equal(as.numeric(sort(precursorMz(aa))[1]),437.80401611) ## [*]
  ## Retention time
  expect_equal(length(rtime(aa)),5)
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_equal(range(rtime(aa)),c(1501.35,1502.31))
  expect_equal(as.numeric(sort(rtime(aa))[1]),1501.35) ## [*]
  ## [*] using as.numeric because rtime and precursorMz return named numerics
  ## Meta data
  expect_equal(dim(fData(aa)),c(5,1))
  expect_equal(dim(pData(aa)),c(0,0))
  ## subsetting
  expect_true(all.equal(aa[["X4"]],assayData(aa)[["X4"]]))
  sub.aa <- aa[1:2]  
  expect_true(all.equal(sub.aa[["X1"]], assayData(sub.aa)[["X1"]]))
  expect_true(all.equal(sub.aa[["X2"]],assayData(sub.aa)[["X2"]]))
  expect_equal(fData(sub.aa),fData(aa)[1:2,,drop=FALSE])
  my.prec <- precursorMz(aa)[1]
  my.prec.aa <- extractPrecSpectra(aa,my.prec)
  expect_true(all(precursorMz(my.prec.aa)==my.prec))
  expect_equal(length(my.prec.aa),2)
  expect_equal(ls(assayData(my.prec.aa)),paste("X",c(1,3),sep=""))
  ## subsetting errors
  expect_error(aa[[1:3]],"subscript out of bounds")
  expect_error(aa[c("X1","X2")],"subsetting works only with numeric or logical")
  expect_error(aa[["AA"]]," object 'AA' not found")
  expect_error(aa[1:10],"subscript out of bounds")
  ## testing that accessors return always attributes in same order
  precMzNames <- names(precursorMz(aa))
  ticNames <- names(tic(aa))
  expect_that(precMzNames,equals(ticNames))
  precChNames <- names(precursorCharge(aa))
  expect_that(precMzNames,equals(precChNames))
  aqnNames <- names(acquisitionNum(aa))
  expect_that(precMzNames,equals(aqnNames))
  rtNames <- names(rtime(aa))
  expect_that(precMzNames,equals(rtNames))
  pkCntNames <- names(peaksCount(aa))
  expect_that(precMzNames,equals(pkCntNames))
  mslNames <- names(msLevel(aa))
  expect_that(precMzNames,equals(mslNames))
  coleNames <- names(collisionEnergy(aa))
  expect_that(precMzNames,equals(coleNames))
  intNames <- names(intensity(aa))
  expect_that(precMzNames,equals(intNames))
  mzNames <- names(mz(aa))
  expect_that(precMzNames,equals(mzNames))
  ffNames <- names(fromFile(aa))
  expect_that(precMzNames,equals(ffNames))  
})

context("MSnExp processing")

test_that("MSnExp processing", {
  ## removeReporters
  expect_true(all.equal(removeReporters(itraqdata[[1]]),itraqdata[[1]]))
  expect_true(all.equal(removeReporters(itraqdata,reporters=iTRAQ4,verbose=FALSE)[[1]],
                        removeReporters(itraqdata[[1]],reporters=iTRAQ4)))
  ## quantitation should be 0
  expect_true(all(quantify(removeReporters(itraqdata[[1]],reporters=iTRAQ4),"max",iTRAQ4)[[1]]==0))
})

## ! Issues with edited dummy data for MS1 uploading, although
## ! things work fine for original data set. Commented these
## ! tests for the moment
## test_that("readMzXMLData and dummy MSnExp msLevel 1 instance", {
##   file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
##   aa <- readMzXMLData(file,msLevel=1,verbose=FALSE)
##   expect_that(class(aa)=="MSnExp",is_true())
##   expect_equal(length(aa),equals(5))
##   ## MS levels
##   expect_that(length(msLevel(aa)),equals(5))
##   expect_that(unique(msLevel(aa)),equals(1))
##   ## Retention time
##   expect_that(length(rtime(aa)),equals(5))
##   expect_that(rtime(aa)[1],is_a("numeric"))
##   expect_that(range(rtime(aa)),
##               equals(c(1982.08,3015.47)))
##   expect_that(as.numeric(polarity(aa)),equals(rep(-1,length(aa)))) ## [*]
##   expect_that(as.numeric(rtime(aa)[1]),equals(1982.08)) ## [*]
##   ## [*] using as.numeric because rtime and precursorMz return named numerics
## })

context("MSnExp data")

test_that("spectra order and integrity", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  clean.aa <- clean(aa,verbose=FALSE)
  rmpeaks.aa <- removePeaks(aa,verbose=FALSE)
  expect_that(ls(assayData(clean.aa)),equals(ls(assayData(aa))))
  expect_that(ls(assayData(rmpeaks.aa)),equals(ls(assayData(aa))))
  int <- c(0,2,3,1,0,0,1)
  sp <- new("Spectrum2",
            intensity = int,
            mz = 1:length(int))
  rsp <- removePeaks(sp)  
  expect_that(peaksCount(sp),equals(length(int)))
  expect_that(tic(sp),equals(sum(int)))
  expect_that(all.equal(removePeaks(sp),rsp),is_true())
  expect_that(tic(removePeaks(sp,1)),equals(6))  
  expect_that(tic(removePeaks(sp,3)),equals(0))
  expect_that(tic(removePeaks(sp,max(intensity(sp)))),equals(0))
  expect_that(peaksCount(sp),equals(peaksCount(rsp)))
  expect_that(peaksCount(clean(rsp)),equals(6))
  expect_that(peaksCount(clean(sp)),equals(7))  
  expect_that(all.equal(removePeaks(sp,0),sp),is_true())
})


test_that("MSnExp normalisation", {
  aa <- itraqdata[1:3]
  bb <- normalise(aa,"max")
  expect_true(all(sapply(intensity(bb),max)==1))
  expect_true(all(unlist(sapply(intensity(aa),order))==unlist(sapply(intensity(bb),order))))
})
