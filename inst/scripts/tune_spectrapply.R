## Purpose of this script:
## - evaluate memory demand by the spectrapply and chromatogram methods.
## - tune spectrapply, eventually removing the gc calls etc.
## - do torture tests to see if changes don't cause/introduce new bugs.
library(MSnbase)
library(microbenchmark)
rm(list = ls())
gc()
gc()

## variables to be protected
protect <- "f"


#####
## Evaluate memory demand:
## - Reading a profile-more mzML file.
f <- "/Users/jo/data/2017/profile/090217_21m_RT_-80_2h_b-AB-profile.mzML.gz"

od <- readMSData2(f)
## macOS 10.12: 489.7 MB

sps <- spectrapply(od)
## macOS 10.12: 1.58 GB
gc()
## macOS 10.12: 963 MB; spectra/spectrapply doesn't seem to release all mem.
rm(list = ls()[!(ls() %in% protect)])
gc()
gc()

im <- readMSData(f, msLevel = 1)
## macOS 10.12: 1.05 GB
gc()
## macOS 10.12: 1.03 GB; gc() doesn't release any additional memory.
rm(list = ls()[!(ls() %in% protect)])
gc()
gc()

## Conclusions from that: spectrapply doesn't seem to free up memory. Actually,
## memory was just not released, i.e. objects were removed but garbage collector
## did just not yet kick in.

## Changes to spectrapply:
rm(list = ls()[!(ls() %in% protect)])
od <- readMSData2(f)

sps <- MSnbase:::spectrapply2(od)
## macOS 10.12: 2.31 GB
gc()
## macOS 10.12: 965 MB


######### BENCHMARKING
microbenchmark(MSnbase:::spectrapply2(od), spectrapply(od), times = 5)
## Unit: seconds
##                        expr      min       lq     mean   median       uq
##  MSnbase:::spectrapply2(od) 24.86483 24.92614 25.03653 25.11650 25.13660
##             spectrapply(od) 25.01711 25.09064 25.32057 25.14135 25.25528
##       max neval cld
##  25.13857     5   a
##  26.09850     5   a
