library(MSnbase)
library(msdata)
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

## readMSData
for (i in 1:10000) {
    if (i %% 200 == 0)
        cat(i, "\n")
    res <- readMSData(f)
    sp <- spectra(res)
    sp <- res[[12]]
}

## readMSData2
for (i in 1:7500) {
    if (i %% 200 == 0)
        cat(i, "\n")
    res <- readMSData2(f, verbose = FALSE)
    sp <- spectra(res)
    sp <- res[[12]]
    res <- readMSData2(mzf, verbose = FALSE)
    sp <- spectra(res)
    sp <- res[[12]]
}
