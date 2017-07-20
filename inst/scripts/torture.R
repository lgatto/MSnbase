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
    res <- readMSData2(f[1], verbose = FALSE)
    sp <- res[[6]]
    res <- readMSData2(f[2], verbose = FALSE)
    sp <- spectra(res)
}

## Issue #170, i.e. random "[MSData::Spectrum::getMZIntensityPairs()] Sizes do
## not match." with the pwiz backend (segfault with ramp). The solution was to
## read call mzR::header() before the mzR::peaks, but this significantly slows
## things down.

library(MSnbase)
SN <- "/Users/jo/data/2016/2016-11/NoSN/"
fl <- dir(SN, full.names = TRUE)
tmp <- readMSData2(fl, msLevel. = 1)

register(SerialParam())
for (i in 1:10) {
    cat("--- ", i, " ---", "\n")
    cat("first spectrapply\n")
    sp <- MSnbase:::spectrapply2(tmp)
    cat("second spectrapply\n")
    sp <- MSnbase:::spectrapply2(tmp)
    tmp <- filterRt(tmp, rt = c(5, 500))
    cat("third spectrapply after filter rt\n")
    sp <- MSnbase:::spectrapply2(tmp)
    cat("\n\n")
}
