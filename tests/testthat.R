# setting R_TESTS to empty string because of
# https://github.com/hadley/testthat/issues/144
# revert this when that issue in R is fixed.
Sys.setenv("R_TESTS" = "")
library("testthat")
library("MSnbase")
setMSnbaseVerbose(FALSE)
## register(SerialParam()) ## see issue 205

## Erwinia
f <- msdata::proteomics(full.names = TRUE,
                        pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
tmt_erwinia_in_mem_ms1 <- readMSData(f, msLevel = 1)
tmt_erwinia_in_mem_ms2 <- readMSData(f, msLevel = 2)
tmt_erwinia_on_disk <- readMSData(f, mode = "onDisk")
tmt_erwinia_on_disk_ms1 <- readMSData(f, msLevel = 1, mode = "onDisk")
tmt_erwinia_on_disk_ms2 <- readMSData(f, msLevel = 2, mode = "onDisk")
## subset by rt
tmt_im_ms1_sub <- filterRt(tmt_erwinia_in_mem_ms1, c(1200, 1250))
tmt_im_ms2_sub <- filterRt(tmt_erwinia_in_mem_ms2, c(1200, 1250))
tmt_od_sub <- filterRt(tmt_erwinia_on_disk, c(1200, 1250))
tmt_od_ms1_sub <- filterRt(tmt_erwinia_on_disk_ms1, c(1200, 1250))
tmt_od_ms2_sub <- filterRt(tmt_erwinia_on_disk_ms2, c(1200, 1250))

## microtofq
f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
       system.file("microtofq/MM8.mzML", package = "msdata"))
microtofq_in_mem_ms1 <- readMSData(f, msLevel = 1)
microtofq_on_disk_ms1 <- readMSData(f, msLevel = 1, mode = "onDisk")
microtofq_on_disk <- readMSData(f, mode = "onDisk")

## extdata mzML
f <- dir(system.file(package = "MSnbase", dir = "extdata"),
         full.name = TRUE, pattern = "mzXML$")
extdata_mzXML_in_mem_ms2 <- readMSData(f, verbose = FALSE, centroided. = FALSE)
extdata_mzXML_on_disk <- readMSData(f, centroided. = FALSE, mode = "onDisk")
extdata_mzXML_on_disk_ms2 <- readMSData(f, msLevel = 2, centroided. = FALSE, mode = "onDisk")

sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
sciex <- readMSData(sf, mode = "onDisk")

sciex_inmem <- readMSnExperiment(sf, backend = BackendMemory())
sciex_mzr <- readMSnExperiment(sf, backend = BackendMzR())
sciex_h5 <- readMSnExperiment(sf, backend = BackendHdf5(),
                              path = file.path(tempdir(), "sciex_h5"))

test_check("MSnbase")
