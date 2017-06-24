# setting R_TESTS to empty string because of
# https://github.com/hadley/testthat/issues/144
# revert this when that issue in R is fixed.
Sys.setenv("R_TESTS" = "")
library("testthat")
library("MSnbase")
setMSnbaseVerbose(FALSE)
register(SerialParam()) ## see issue 205

## Erwinia
f <- msdata::proteomics(full.names = TRUE,
                        pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
tmt_erwinia_in_mem_ms1 <- readMSData(f, msLevel = 1)
tmt_erwinia_in_mem_ms2 <- readMSData(f, msLevel = 2)
tmt_erwinia_on_disk <- readMSData2(f)
tmt_erwinia_on_disk_ms1 <- readMSData2(f, msLevel = 1)
tmt_erwinia_on_disk_ms2 <- readMSData2(f, msLevel = 2)

## microtofq
f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
       system.file("microtofq/MM8.mzML", package = "msdata"))
microtofq_in_mem_ms1 <- readMSData(f, msLevel = 1)
microtofq_on_disk_ms1 <- readMSData2(f, msLevel = 1)
microtofq_on_disk <- readMSData2(f)

## extdata mzML
f <- dir(system.file(package = "MSnbase", dir = "extdata"),
         full.name = TRUE, pattern = "mzXML$")
extdata_mzXML_in_mem_ms2 <- readMSData(f, verbose = FALSE, centroided. = FALSE)
extdata_mzXML_on_disk <- readMSData2(f, centroided. = FALSE)
extdata_mzXML_on_disk_ms2 <- readMSData2(f, msLevel = 2, centroided. = FALSE)


test_check("MSnbase")
