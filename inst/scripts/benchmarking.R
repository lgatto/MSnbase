## Load libs
library(MSnbase)
library(testthat)
library(microbenchmark)

## Define test files
library(msdata)
fms1 <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
          system.file("microtofq/MM8.mzML", package = "msdata"))
fmsn <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

## Read files.
multiFilesInMem <- readMSData(fms1, msLevel. = 1, centroided. = TRUE)
multiFilesOnDisk <- readMSData(fms1, centroided. = TRUE, mode = "onDisk")

multiMsInMem1 <- readMSData(fmsn, msLevel. = 1)
multiMsInMem2 <- readMSData(fmsn, msLevel. = 2)
multiMsOnDisk <- readMSData(fmsn, mode = "onDisk")
multiMsOnDisk1 <- filterMsLevel(multiMsOnDisk, msLevel. = 1)

## Extract all spectra.
microbenchmark(spectra(multiFilesInMem), spectra(multiFilesOnDisk), times = 10)
microbenchmark(spectra(multiMsInMem1),
               spectra(filterMsLevel(multiMsOnDisk, msLevel. = 1)), times = 10)

## filterFile
microbenchmark(filterFile(multiFilesInMem, file = 2),
               filterFile(multiFilesOnDisk, file = 2), times = 10)

## Extract single spectrum.
microbenchmark(multiMsInMem1[[23]], multiMsOnDisk1[[23]], times = 10)

## Subset data.
microbenchmark(multiMsInMem1[7:23], multiMsOnDisk1[7:23], times = 10)


## Setting centroided to TRUE
microbenchmark(centroided(multiMsInMem1) <- TRUE,
               centroided(multiMsOnDisk1) <- TRUE, times = 10)
