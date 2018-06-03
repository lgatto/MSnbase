library(MSnbase)
library(msdata)
f <- msdata::proteomics(full.names = TRUE)

## readMSData
for (i in 1:10000) {
    if (i %% 200 == 0)
        cat(i, "\n")
    res <- readMSData(f[1], verbose = FALSE, msLevel. = 1)
    res <- readMSData(f[1], verbose = FALSE, msLevel. = 2)
    sp <- res[[12]]
    sp <- spectra(res)
}

## readMSData onDisk
for (i in 1:7500) {
    if (i %% 200 == 0)
        cat(i, "\n")
    res <- readMSData(f[4], verbose = FALSE, mode = "onDisk")
    sp <- res[[6]]
    res <- readMSData(f[4], verbose = FALSE, mode = "onDisk")
    sp <- spectra(res)
    chr <- chromatogram(res)
}
## OK without gc()

## spectrapply2
res <- readMSData(f[4], verbose = FALSE, mode = "onDisk")
for (i in 1:10000) {
    if (i %% 200 == 0)
        cat(i, "\n")
    ## sp <- res[[6]]
    sp <- spectrapply(res)
    sp <- spectrapply(res, FUN = mz)
    tmp <- filterRt(res, rt = c(1108, 1109))
    sp <- spectrapply(tmp)
}
## OK without gc()

## spectrapply2
res <- readMSData(f[4], verbose = FALSE, mode = "onDisk")
for (i in 1:10000) {
    if (i %% 200 == 0)
        cat(i, "\n")
    ## sp <- res[[6]]
    res <- readMSData(f[4], verbose = FALSE, mode = "onDisk")
    sp <- MSnbase:::spectrapply2(res)
    sp <- MSnbase:::spectrapply2(res, FUN = mz)
    tmp <- filterRt(res, rt = c(1108, 1109))
    sp <- MSnbase:::spectrapply2(tmp)
}
## OK without gc()

## Issue #170, i.e. random "[MSData::Spectrum::getMZIntensityPairs()] Sizes do
## not match." with the pwiz backend (segfault with ramp). The solution was to
## read call mzR::header() before the mzR::peaks, but this significantly slows
## things down.

library(MSnbase)

torturing <- function(x) {
    tmp <- readMSData(x, msLevel. = 1, mode = "onDisk")
    register(SerialParam())
    for (i in 1:10) {
        cat("--- ", i, " ---", "\n")
        cat("first spectrapply\n")
        sp <- MSnbase::spectrapply(tmp, FUN = function(z) {max(mz(z))})
        ## sp has a size of 21.9GB!!!
        rm(sp)
        gc()
        cat("second spectrapply\n")
        sp <- MSnbase::spectrapply(tmp, FUN = function(z) {max(mz(z))})
        rm(sp)
        gc()
        tmp <- filterRt(tmp, rt = c(5, 500))
        cat("third spectrapply after filter rt\n")
        sp <- MSnbase::spectrapply(tmp, FUN = function(z) {max(mz(z))})
        cat("\n\n")
    }
}

torturing2 <- function(x) {
    tmp <- readMSData(x, msLevel. = 1, mode = "onDisk")
    register(SerialParam())
    for (i in 1:10) {
        cat("--- ", i, " ---", "\n")
        cat("first spectrapply\n")
        sp <- MSnbase:::spectrapply2(tmp, FUN = function(z) {max(mz(z))})
        ## sp has a size of 21.9GB!!!
        rm(sp)
        gc()
        cat("second spectrapply\n")
        sp <- MSnbase:::spectrapply2(tmp, FUN = function(z) {max(mz(z))})
        rm(sp)
        gc()
        tmp <- filterRt(tmp, rt = c(5, 500))
        cat("third spectrapply after filter rt\n")
        sp <- MSnbase:::spectrapply2(tmp, FUN = function(z) {max(mz(z))})
        cat("\n\n")
    }
}

## All tests are now WITHOUT gc().

SN <- "/Users/jo/data/2016/2016-11/NoSN/"
fl <- dir(SN, full.names = TRUE)

torturing(fl)
## macOS: 2x FAIL
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 2x OK.

## Windows: 2x OK.

setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS without gc(): 1x OK

## spectrapply2
setMSnbaseFastLoad(FALSE)
torturing2(fl)
## macOS without gc(): 1x OK

fl <- dir("/Users/jo/data/2017/2017_02/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing(fl)
## macOS: 2x FAIL
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 2x OK.

## Windows: 2x OK.


fl <- dir("/Users/jo/data/2017/2017_02/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 2x OK.

fl <- dir("/Users/jo/data/2017/2017_02/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing2(fl)
## maxOS: 2x OK.


fl <- dir("/Users/jo/data/2017/nalden01/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing(fl)
## macOS: 2x FAIL
## ---  8  --- 
## first spectrapply
## second spectrapply
## third spectrapply after filter rt
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.
## ---  6  --- 
## first spectrapply
## second spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 2x OK.
## Windows: 2x OK.

fl <- dir("/Users/jo/data/2017/nalden01/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 1x FAIL: that's problematic; it's using the header(last spectrum).
## ---  3  --- 
## first spectrapply
## second spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## macOS without gc(): 1x FAIL, like above!
## Instead of reading the last header we read ALL header for spectrapply
## macOS with header(all): 1x OK, 1x FAILED (like above).

fl <- dir("/Users/jo/data/2017/nalden01/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing2(fl)
## macOS without gc(): 3x OK.


fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing(fl)
## macOS: 1x FAIL:
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 2x OK.

## Windows: 2x OK.

fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 2x OK.

fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing2(fl)
## macOS: 1x FAIL (typical error).

fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing2(fl)
## macOS: 2x OK



###################
## Testing eventual less memory demanding versions of spectrapply.
## Variangs are:
## o spectrapply: reading all mz/intensity values. Creating the list of
##   Spectrum objects for them and calling lapply on the function.
## o spectrapply3: processing each spectrum separately, i.e. load the data,
##   create the spectrum object and apply the function.
## o spectrapply4: same as _3, but using a for loop instead of split, lapply

## What I should test here:
## 1) Can I use fastLoad = TRUE for spectrapply4? NO!
## 2) Do I need the gc() for the spectrapply4? TEST LATER.

library(MSnbase)
register(SerialParam())
setMSnbaseFastLoad(FALSE)

fl <- dir("/Users/jo/data/2016/2016-11/NoSN", full.names = TRUE)
tmp <- readMSData(fl, mode = "onDisk")

for (i in 1:10) {
    cat("\nIteration", i, "\n\n")
    maxi <- MSnbase:::spectrapply2(tmp, FUN = function(z) {max(mz(z))})
}

## 3) Compare the performance.LLLLLLL
library(MSnbase)
library(microbenchmark)
library(testthat)

## mzML
fl <- system.file("microtofq/MM14.mzML", package = "msdata")
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
sp_1 <- spectrapply(od, FUN = mz)
sp_2 <- MSnbase:::spectrapply2(od, FUN = mz)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr      min       lq      mean   median        uq
##             spectrapply(od) 85.15372 88.22235 109.80866 88.24776  89.24869
##  MSnbase:::spectrapply2(od) 95.21179 96.86848  99.96095 98.65392 103.79570
##       max neval cld
##  198.1708     5   a
##  105.2749     5   a

## With a method
microbenchmark(spectrapply(od, FUN = mz), MSnbase:::spectrapply2(od, FUN = mz),
               times = 5)
## Unit: milliseconds
##                                  expr      min        lq      mean    median
##             spectrapply(od, FUN = mz) 82.87224  95.75993  96.40574  97.41207
##  MSnbase:::spectrapply2(od, FUN = mz) 96.75382 101.77152 104.28466 101.86804
##         uq      max neval cld
##   99.92526 106.0592     5   a
##  109.38800 111.6419     5   a

## Large file:
fl <- "/Users/jo/data/2017/2017_02/090217_111m_RT_-80_24h_a.mzML"
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
setMSnbaseFastLoad(FALSE)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr       min        lq     mean   median       uq
##             spectrapply(od)  969.9947  992.4444 1161.431 1166.618 1187.340
##  MSnbase:::spectrapply2(od) 1177.2195 1190.1962 1345.919 1380.541 1458.512
##       max neval cld
##  1490.758     5   a
##  1523.124     5   a

microbenchmark(spectrapply(od, FUN = mz), MSnbase:::spectrapply2(od, FUN = mz),
               times = 5)
## Unit: milliseconds
##                                  expr       min        lq     mean   median
##             spectrapply(od, FUN = mz)  970.0606  986.7927 1033.393 1015.398
##  MSnbase:::spectrapply2(od, FUN = mz) 1163.2397 1181.7938 1291.004 1197.522
##        uq      max neval cld
##  1017.076 1177.636     5  a 
##  1384.279 1528.184     5   b

## gzipped mzML
fl <- system.file("proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML.gz", package = "msdata")
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 23.77440 24.01636 25.46504 24.95036 27.19821
##  MSnbase:::spectrapply2(od) 20.35279 20.77263 21.26155 21.26453 21.90014
##       max neval cld
##  27.38587     5   b
##  22.01765     5  a 


## mzXML
fl <- system.file("lockmass/LockMass_test.mzXML", package = "msdata")
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 177.9468 181.9812 195.6146 186.8247 205.5001
##  MSnbase:::spectrapply2(od) 224.5103 259.7205 279.5995 271.8644 313.5613
##       max neval cld
##  225.8203     5  a 
##  328.3411     5   b

## At last with the same file but gzipped...
fl <- "/Users/jo/data/2017/mzXML/1405_blk1.mzXML.gz"
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 41.49147 41.94676 43.41425 42.21200 42.23053
##  MSnbase:::spectrapply2(od) 39.08837 43.56239 44.41015 45.10258 47.06231
##       max neval cld
##  49.19049     5   a
##  47.23510     5   a

fl <- system.file("cdf/ko15.CDF", package = "msdata")
od <- readMSData(fl, mode = "onDisk")
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr        min         lq       mean     median
##             spectrapply(od)   266.9588   268.5773   338.2803   269.8898
##  MSnbase:::spectrapply2(od) 13483.9753 13582.1238 13843.0532 13800.6454
##         uq        max neval cld
##    402.648   483.3276     5  a 
##  14100.982 14247.5391     5   b
