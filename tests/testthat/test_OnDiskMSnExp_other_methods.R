context("OnDiskMSnExp class, other methods")

############################################################
## Load the required data files.
.getMzMLFiles <- function() {
    ## Return the mzML files, the ones from the XXX package, or if run
    ## locally, some of my test files.
    HOST <- unlist(strsplit(system("hostname", intern = TRUE), split = ".",
                            perl = FALSE, fixed = TRUE))[1]
    if (HOST == "macbookjo") {
        mzfiles <- dir("/Users/jo/R-workspaces/EURAC/2016/2016-04-21-PolarMetabolom/data/mzML/",
                       pattern = "POS_C_O", full.names = TRUE)
    } else {
        require(msdata)
        mzfiles <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                     system.file("microtofq/MM8.mzML", package = "msdata"))
    }
    return(mzfiles)
}
mzf <- .getMzMLFiles()[1:2]
## Load the data as an MSnExp into memory.
inMem <- readMSData(files = mzf, msLevel = 1, centroided = TRUE)
## Load the data as OnDiskMSnExp.
suppressWarnings(
    onDisk <- readMSData2(files = mzf, msLevel = 1, centroided = TRUE)
)

## Read another mzML file.
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem2 <- readMSData(f, centroided = NA, verbose = FALSE)  ## That's the MS 2 data.
inmem1 <- readMSData(f, centroided = NA, verbose = FALSE, msLevel = 1)  ## MS 1 data.
ondisk <- readMSData2(f, verbose = FALSE)
ondisk1 <- readMSData2(f, msLevel = 1, verbose = FALSE)
ondisk2 <- readMSData2(f, msLevel = 2, verbose = FALSE)


############################################################
## plot
test_that("OnDiskMSnExp plot", {
    ## Would be nice to know what the plot function is actually doing though...
    ## seems I can forget that for larger experiments; takes way to long.
})

############################################################
## trimMz
test_that("OnDiskMSnExp trimMz", {
    ## Comparing timings and results for the trimMz.
    inMem <- inmem1
    onDisk <- ondisk1
    system.time(
        inMemMz <- trimMz(inMem, mzlim=c(500, 550))
    ) ## 7.3 sec.
    system.time(
        onDiskMz <- trimMz(onDisk, mzlim=c(500, 550))
    ) ## woah, 0.009 sec (what a surprise ;) )
    expect_true(all.equal(inMemMz, onDiskMz))
})

############################################################
## normalize
test_that("Compare OnDiskMSnExp and MSnExp normalize", {
    ## Compare timings and results for normalize.
    system.time(
        inMemN <- normalize(inMem)
    )  ## 14.2 sec
    system.time(
        onDiskN <- normalize(onDisk)
    )  ## 0.005
    ## Get and compare spectra.
    expect_true(all.equal(inMemN, onDiskN))

    ## Compare intensity values.
    expect_identical(intensity(inMemN), intensity(onDiskN))
    ## Compare sub-setted result.
    inMemNSub <- inMemN[1:20]
    onDiskNSub <- onDiskN[1:20]
    expect_true(all.equal(inMemNSub, onDiskNSub))
})

############################################################
## smooth
test_that("Compare OnDiskMSnExp and MSnExp smooth", {

    ## The same for MSn data
    ondisk <- smooth(ondisk)
    suppressWarnings(
        inmem1 <- smooth(inmem1)
    )
    suppressWarnings(
        inmem2 <- smooth(inmem2)
    )
    suppressWarnings(
        expect_true(all.equal(filterMsLevel(ondisk, msLevel. = 1), inmem1))
    )
    suppressWarnings(
        expect_true(all.equal(filterMsLevel(ondisk, msLevel. = 2), inmem2))
    )
})

############################################################
## compareSpectra
test_that("Compare OnDiskMSnExp and MSnExp compareSpectra", {
    system.time(
        csp <- compareSpectra(inmem1)
    ) ## 39.8 sec
    system.time(
        csp2 <- compareSpectra(ondisk1)
    )
    rownames(csp) <- colnames(csp) <- NULL
    rownames(csp2) <- colnames(csp2) <- NULL
    expect_identical(csp, csp2)
})

