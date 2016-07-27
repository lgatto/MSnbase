context("OnDiskMSnExp class, other methods")

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
mzf <- .getMzMLFiles(TRUE)[1:2]
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

############################################################
## pickPeaks
test_that("Compare OnDiskMSnExp and MSnExp pickPeaks", {
    ## Setting centroided FALSE to avoid warnings...
    centroided(inmem1) <- FALSE
    system.time(
        pp <- pickPeaks(inmem1)
    ) ## 0.5 sec
    centroided(ondisk1) <- FALSE
    pp2 <- pickPeaks(ondisk1)
    expect_true(all.equal(pp, pp2))
    pc <- peaksCount(pp)
    pc2 <- peaksCount(pp2)
    names(pc) <- NULL
    names(pc2) <- NULL
    expect_identical(pc, pc2)
})

############################################################
## estimateNoise
test_that("Compare OnDiskMSnExp and MSnExp estimateNoise", {
    centroided(inmem1) <- FALSE
    system.time(
        en <- estimateNoise(inmem1)
    ) ## 0.142
    centroided(ondisk1) <- FALSE
    system.time(
        en2 <- estimateNoise(ondisk1)
    ) ## 1.1
    names(en) <- NULL
    names(en2) <- NULL
    expect_identical(en, en2)
    ## Same with method = SuperSmoother
    centroided(inmem1) <- FALSE
    system.time(
        en <- estimateNoise(inmem1, method = "SuperSmoother")
    ) ## 1.149
    centroided(ondisk1) <- FALSE
    system.time(
        en2 <- estimateNoise(ondisk1, method = "SuperSmoother")
    ) ## 2.9
    names(en) <- NULL
    names(en2) <- NULL
    expect_identical(en, en2)
})

############################################################
## extractPrecSpectra
test_that("Compare OnDiskMSnExp and MSnExp extractPrecSpectra", {
    precMzs <- precursorMz(inmem2)
    exP <- sample(1:length(precMzs), 20)
    extsp <- extractPrecSpectra(inmem2, prec = precMzs[exP])
    extsp2 <- extractPrecSpectra(ondisk2, prec = precMzs[exP])

    expect_true(all.equal(extsp, extsp2))

    ## Now with a multi-MS-level object.
    extsp3 <- extractPrecSpectra(ondisk, prec = precMzs[exP])
    expect_true(all.equal(extsp2, extsp3))

    ## Check that the precursorMz matches:
    expect_identical(unname(precMzs[sort(exP)]),
                     unname(precursorMz(extsp2)))
})

############################################################
## isCentroided
test_that("isCentroided on OnDiskMSnExp", {

    expect_true(all(isCentroided(onDisk, verbose = FALSE)))
    expect_true(all(isCentroided(ondisk2, verbose = FALSE)))

})

############################################################
## Test precursor* methods
test_that("Test precursor* for OnDiskMSnExp", {
    ## o precursorMz
    pmz <- precursorMz(ondisk)
    pmz2 <- precursorMz(ondisk2)
    expect_true(all(is.na(pmz[msLevel(ondisk) == 1])))
    expect_identical(pmz2, pmz[names(pmz2)])

    ##  Finally compare to inmem.
    pmz <- precursorMz(inmem2)
    names(pmz) <- names(pmz2) <- NULL
    expect_identical(pmz, pmz2)

    ## o precursorCharge
    pch <- precursorCharge(ondisk)
    pch2 <- precursorCharge(ondisk2)
    expect_true(all(is.na(pch[msLevel(ondisk) == 1])))
    expect_identical(pch2, pch[names(pch2)])

    ##  Finally compare to inmem.
    pch <- precursorCharge(inmem2)
    names(pch) <- names(pch2) <- NULL
    expect_identical(pch, pch2)

    ## o precursorIntensity
    pint <- precursorIntensity(ondisk)
    pint2 <- precursorIntensity(ondisk2)
    expect_true(all(is.na(pint[msLevel(ondisk) == 1])))
    expect_identical(pint2, pint[names(pint2)])

    ##  Finally compare to inmem.
    pint <- precursorIntensity(inmem2)
    names(pint) <- names(pint2) <- NULL
    expect_identical(pint, pint2)

    ## o precScanNum
    pcn <- precScanNum(ondisk)
    pcn2 <- precScanNum(ondisk2)
    expect_true(all(is.na(pcn[msLevel(ondisk) == 1])))
    expect_identical(pcn2, pcn[names(pcn2)])

    ##  Finally compare to inmem.
    pcn <- precScanNum(inmem2)
    names(pcn) <- names(pcn2) <- NULL
    expect_identical(pcn, pcn2)
})

