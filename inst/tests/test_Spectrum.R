context("Spectrum class")

test_that("Spectrum validity", {
  expect_true(validObject(new("Spectrum1")))
  expect_true(validObject(new("Spectrum2")))
})


test_that("Spectrum processing", {
  int <- c(0,1,2,3,1,0,0,0,0,1,3,10,6,2,1,0,1,2,0,0,1,5,10,5,1)
  sp <- new("Spectrum2",
            intensity=int,
            mz=1:length(int))
  ## removePeaks
  defaultT <- min(intensity(sp)[intensity(sp)>0])
  sp2a <- removePeaks(sp)
  sp2b <- removePeaks(sp,defaultT)
  sp2c <- removePeaks(sp,3)
  expect_that(identical(sp2a,sp2b),is_true())
  expect_that(identical(sp,sp2b),is_true())
  expect_that(identical(sp,sp2c),is_false())
  expect_that(tic(sp),equals(sum(int)))
  expect_that(peaksCount(sp2c),equals(peaksCount(sp)))
  expect_that(tic(sp),equals(55))
  expect_that(tic(sp2c),equals(45))
  ## clean
  sp3 <- clean(sp)
  expect_that(tic(sp),equals(tic(sp3)))
  expect_that(peaksCount(sp),equals(length(int)))
  expect_that(peaksCount(sp3),equals(23))
  ##trimMz
  sp4 <- trimMz(sp,c(10,20))
  expect_that(intensity(sp4),equals(int[10:20]))
  expect_that(mz(sp4),equals(10:20))
  expect_that(peaksCount(sp4),equals(length(10:20)))
  expect_that(tic(sp4),equals(sum(int[10:20])))  
})


test_that("Spectrum quantification", {
  ## dummy Spectrum
  int <- c(0,2,3,1,0)
  mz <- c(114.11,
          114.12,
          114.13,
          114.14,
          114.145)
  sp <- new("Spectrum2",
            intensity=int,
            mz=mz)
  expect_true(validObject(sp))
  expect_equal(MSnbase:::getCurveWidth(sp,iTRAQ4[1]),list(lwr=1,upr=5))
  expect_equal(as.numeric(quantify(sp,"sum",iTRAQ4[1])$peakQuant),6)
  expect_equal(as.numeric(quantify(sp,"max",iTRAQ4[1])$peakQuant),3)
  expect_that(as.numeric(quantify(sp,"trap",iTRAQ4[1])$peakQuant),
              equals((0.01*2)/2+
                     (0.01*2)  +
                     (0.01*1)/2+
                     0.01*1    +
                     (0.01*2)/2+
                     (0.01*0.5)/2))
  print("Warnings expected because there is not data for iTRAQ4[2].")
  expect_true(as.logical(is.na(quantify(sp,"sum",iTRAQ4[2])$peakQuant)))
  expect_warning(quantify(sp,"sum",iTRAQ4[2])$peakQuant)
})

