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

## readMSData2
for (i in 1:7500) {
    if (i %% 200 == 0)
        cat(i, "\n")
    res <- readMSData2(f[4], verbose = FALSE)
    sp <- res[[6]]
    res <- readMSData2(f[4], verbose = FALSE)
    sp <- spectra(res)
}
## OK without gc()

## spectrapply2
res <- readMSData2(f[4], verbose = FALSE)
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
res <- readMSData2(f[4], verbose = FALSE)
for (i in 1:10000) {
    if (i %% 200 == 0)
        cat(i, "\n")
    ## sp <- res[[6]]
    res <- readMSData2(f[4], verbose = FALSE)
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
    tmp <- readMSData2(x, msLevel. = 1)
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

SN <- "/Users/jo/data/2016/2016-11/NoSN/"
fl <- dir(SN, full.names = TRUE)

torturing(fl)
## macOS: 2x FAIL
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 1x FAIL:
## ---  1  ---
## first spectrapply
## Killed

setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 1x FAIL
## ---  1  --- 
## first spectrapply
## Killed: 9

## Linux: 1x FAIL
## ---  1  --- 
## first spectrapply
## Killed: 9

## macOS without gc(): 1x RUNNING

fl <- dir("/Users/jo/data/2017/2017_02/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing(fl)
## macOS: 2x FAIL
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## Linux: 1x FAIL:
## ---  1  ---
## first spectrapply
## second spectrapply
## Killed


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

fl <- dir("/Users/jo/data/2017/nalden01/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 1x FAIL: that's problematic; it's using the header(last spectrum).
## ---  3  --- 
## first spectrapply
## second spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

## macOS without gc(): 1x RUNNING

fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(TRUE)
torturing(fl)
## macOS: 1x FAIL:
## ---  1  --- 
## first spectrapply
## Error in object@backend$getPeakList(x) : 
##   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match.

fl <- dir("/Users/jo/data/2016/2016_06/", full.names = TRUE)
setMSnbaseFastLoad(FALSE)
torturing(fl)
## macOS: 1x FAIL
## ---  1  --- 
## first spectrapply
## second spectrapply
## Killed: 9



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
tmp <- readMSData2(fl)

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
od <- readMSData2(fl)
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
sp_1 <- spectrapply(od, FUN = mz)
sp_2 <- MSnbase:::spectrapply2(od, FUN = mz)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 166.3859 196.2183 191.8018 196.4440 198.3084
##  MSnbase:::spectrapply2(od) 204.4279 206.6654 207.9754 208.4293 209.3346
##       max neval cld
##  201.6526     5  a 
##  211.0201     5   b

## With a method
microbenchmark(spectrapply(od, FUN = mz), MSnbase:::spectrapply2(od, FUN = mz),
               times = 5)
## Unit: milliseconds
##                                  expr      min       lq     mean   median
##             spectrapply(od, FUN = mz) 259.0931 259.4859 260.4230 259.5214
##  MSnbase:::spectrapply2(od, FUN = mz) 269.3629 273.4123 276.2432 275.1020
##        uq      max neval cld
##  260.4673 263.5475     5  a 
##  277.9313 285.4074     5   b

## Large file:
fl <- "/Users/jo/data/2017/2017_02/090217_111m_RT_-80_24h_a.mzML"
od <- readMSData2(fl)
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
setMSnbaseFastLoad(FALSE)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 1.186034 1.193843 1.289523 1.205818 1.411478
##  MSnbase:::spectrapply2(od) 1.337296 1.350252 1.369661 1.371066 1.378367
##       max neval cld
##  1.450443     5   a
##  1.411325     5   a

microbenchmark(spectrapply(od, FUN = mz), MSnbase:::spectrapply2(od, FUN = mz),
               times = 5)
## Unit: seconds
##                                  expr      min       lq     mean   median
##             spectrapply(od, FUN = mz) 1.184725 1.208542 1.207812 1.208744
##  MSnbase:::spectrapply2(od, FUN = mz) 1.333628 1.355051 1.380853 1.391359
##        uq      max neval cld
##  1.217936 1.219111     5  a 
##  1.398276 1.425948     5   b

## gzipped mzML
fl <- system.file("proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML.gz", package = "msdata")
od <- readMSData2(fl)
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 22.16420 22.52596 22.55925 22.55297 22.75131
##  MSnbase:::spectrapply2(od) 19.08312 19.69682 19.94403 20.26291 20.29765
##       max neval cld
##  22.80183     5   b
##  20.37965     5  a 

## mzXML
fl <- system.file("lockmass/LockMass_test.mzXML", package = "msdata")
od <- readMSData2(fl)
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: milliseconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 326.3817 334.5355 341.2993 346.4015 349.5263
##  MSnbase:::spectrapply2(od) 377.3414 383.2008 385.2446 387.4712 388.8911
##       max neval cld
##  349.6515     5  a 
##  389.3184     5   b

## At last with the same file but gzipped...
fl <- "/Users/jo/data/2017/mzXML/1405_blk1.mzXML.gz"
od <- readMSData2(fl)
sp_1 <- spectrapply(od)
sp_2 <- MSnbase:::spectrapply2(od)
expect_equal(sp_1, sp_2)
microbenchmark(spectrapply(od), MSnbase:::spectrapply2(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##             spectrapply(od) 28.92022 28.94184 29.23785 29.20040 29.38479
##  MSnbase:::spectrapply2(od) 29.15081 29.16656 29.43833 29.17894 29.56303
##       max neval cld
##  29.74202     5   a
##  30.13233     5   a
