context("OnDiskMSnExp benchmarks")

## These functions are for benchmarking, not for testing/evaluation of functionality.
############################################################
## Load the required data files.
.getMzMLFiles <- function(force.msdata = FALSE) {
    ## Return the mzML files, the ones from the XXX package, or if run
    ## locally, some of my test files.
    HOST <- unlist(strsplit(system("hostname", intern = TRUE), split = ".",
                            perl = FALSE, fixed = TRUE))[1]
    if (HOST == "macbookjo" & !force.msdata) {
        mzfiles <- dir("/Users/jo/R-workspaces/EURAC/2016/2016-04-21-PolarMetabolom/data/mzML/",
                       pattern = "POS_C_O", full.names = TRUE)
    } else {
        require(msdata)
        mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                     system.file("microtofq/MM8.mzML", package = "msdata"))
    }
    return(mzfiles)
}

mzf <- .getMzMLFiles(force.msdata = TRUE)[1:2]
suppressWarnings(
    odmse <- readMSData2(files = mzf, centroided = TRUE)
)
## Simple benchmark to evaluate whether parallel processing brings a performance gain.
.bench_intensity_serial_parallel <- function() {
    ## The two files have different number of spectra etc, thus we're extracting the
    ## first 100.
    ## For small files (have about 140 spectra per file):
    mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                 system.file("microtofq/MM8.mzML", package = "msdata"))
    suppressWarnings(
        odmse <- readMSData2(files = mzfiles, centroided = TRUE)
    )
    library(microbenchmark)
    microbenchmark(intensity(odmse[1:100], BPPARAM = SerialParam()),
                   intensity(odmse[1:100], BPPARAM = bpparam()), times = 20)

    ## Check if we get what we expect by comparing it to the MSnExp:
    mse <- readMSData(files = mzfiles, centroided = TRUE, msLevel = 1, backend = "ram")
    ints_1 <- intensity(mse[1:100])
    ints_2 <- intensity(odmse[1:100], BPPARAM = SerialParam())
    ints_3 <- intensity(odmse[1:100], BPPARAM = bpparam())
    expect_identical(ints_1, ints_2)
    expect_identical(ints_1, ints_3)

    ## And with large files (3670 spectra per file)
    mzfiles <- .getMzMLFiles()[1:2]
    suppressWarnings(
        odmse <- readMSData2(files = mzfiles, centroided = TRUE)
    )
    ## Get 1000
    microbenchmark(intensity(odmse[1:1000], BPPARAM = SerialParam()),
                   intensity(odmse[1:1000], BPPARAM = bpparam()),
                   times = 20)

    ## Get 2000
    microbenchmark(intensity(odmse[1:2000], BPPARAM = SerialParam()),
                   intensity(odmse[1:2000], BPPARAM = bpparam()),
                   times = 20)

    ## Get 3000
}

