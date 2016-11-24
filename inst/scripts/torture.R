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

