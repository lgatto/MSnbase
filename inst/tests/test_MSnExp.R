context("MSnExp class")

test_that("readMzXMLData and dummy MSnExp msLevel 2 instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  expect_that(class(aa)=="MSnExp",is_true())
  ## checking slots and methods
  expect_that(length(aa),equals(70))
  expect_that(nrow(header(aa)),equals(length(aa)))
  expect_that(names(header(aa)),
              equals(c("index","file","retention.time",
                       "precursor.mz","peaks.count","tic",
                       "ms.level","charge","collision.energy")))
  ## MS levels
  expect_that(length(msLevel(aa)),equals(70))
  expect_that(unique(msLevel(aa)),equals(2))
  expect_that(length(ms1scan(aa)),equals(70))
  expect_that(length(unique(ms1scan(aa))),equals(5))
  ## Precursor MZ
  expect_that(length(precursorMz(aa)),equals(70))
  expect_that(precursorMz(aa)[1],is_a("numeric"))
  expect_that(length(unique(precursorMz(aa))),equals(11))
  expect_that(range(precursorMz(aa)),
              equals(c(424.76651001,630.33178711)))
  expect_that(as.numeric(sort(precursorMz(aa))[1]),equals(424.76651001)) ## [*]
  ## Retention time
  expect_that(length(rtime(aa)),equals(70))
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_that(range(rtime(aa)),
              equals(c(1982.92,3018.43)))
  expect_that(as.numeric(sort(rtime(aa))[1]),equals(1982.92)) ## [*]
  ## [*] using as.numeric because rtime and precursorMz return named numerics
  ## Meta data
  expect_that(dim(fData(aa)),equals(c(70,1)))
  expect_that(dim(pData(aa)),equals(c(0,0)))
  ## subsetting
  expect_that(all.equal(aa[["X64"]],assayData(aa)[["X64"]]),
              is_true())
  sub.aa <- aa[1:2]  
  expect_that(all.equal(sub.aa[["X1"]], assayData(sub.aa)[["X1"]]),is_true())
  expect_that(all.equal(sub.aa[["X10"]],assayData(sub.aa)[["X10"]]),is_true())
  expect_that(fData(sub.aa), equals(fData(aa)[1:2,,drop=FALSE]))
  my.prec <- precursorMz(aa)[1]
  my.prec.aa <- extractPrecSpectra(aa,my.prec)
  expect_that(all(precursorMz(my.prec.aa)==my.prec),is_true())
  expect_that(length(my.prec.aa),equals(4))
  expect_that(ls(assayData(my.prec.aa)),
              equals(c("X65","X67","X69","X70")))
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

test_that("readMzXMLData and dummy MSnExp msLevel 1 instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,msLevel=1,verbose=FALSE)
  expect_that(class(aa)=="MSnExp",is_true())
  expect_that(length(aa),equals(5))
  ## MS levels
  expect_that(length(msLevel(aa)),equals(5))
  expect_that(unique(msLevel(aa)),equals(1))
  ## Retention time
  expect_that(length(rtime(aa)),equals(5))
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_that(range(rtime(aa)),
              equals(c(1982.08,3015.47)))
  expect_that(as.numeric(polarity(aa)),equals(rep(-1,length(aa)))) ## [*]
  expect_that(as.numeric(rtime(aa)[1]),equals(1982.08)) ## [*]
  ## [*] using as.numeric because rtime and precursorMz return named numerics
})


context("data integrity")

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
  expect_that(peaksCount(clean(rsp)),equals(5))
  expect_that(peaksCount(clean(sp)),equals(6))  
  expect_that(all.equal(removePeaks(sp,0),sp),is_true())
})


context("quantification and MSnSet instance")

test_that("quantification", {
  ## dummy Spectrum
  int <- c(0,2,3,1,0)
  mz <- c(114.11,
          114.12,
          114.13,
          114.14,
          114.15)
  sp <- new("Spectrum2",
            intensity=int,
            mz=mz)
  data(iTRAQ4)
  expect_that(validObject(sp),is_true())
  expect_that(getCurveWidth(sp,iTRAQ4[1]),
              equals(list(lwr=1,upr=5)))
  expect_that(as.numeric(quantify(sp,iTRAQ4[1],"sum")$peakQuant),
              equals(6))
  expect_that(as.numeric(quantify(sp,iTRAQ4[1],"max")$peakQuant),
              equals(3))
  expect_that(as.numeric(quantify(sp,iTRAQ4[1],"trap")$peakQuant),
              equals((0.01*2)/2+
                     (0.01*2)  +
                     (0.01*1)/2+
                     0.01*1    +
                     (0.01*2)/2+
                     (0.01*1)/2))
})
