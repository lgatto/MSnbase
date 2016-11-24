context("OnDiskMSnExp class, other methods")

library("msdata")
mzf <-  c(system.file("microtofq/MM14.mzML", package = "msdata"),
          system.file("microtofq/MM8.mzML", package = "msdata"))

inMem <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE,
                    verbose = FALSE)
onDisk <- readMSData2(files = mzf, msLevel. = 1, centroided. = TRUE,
                      verbose = FALSE)

f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia_1")
inmem2 <- readMSData(f, centroided. = NA, verbose = FALSE)  ## That's the MS 2 data.
inmem1 <- readMSData(f, centroided. = NA, verbose = FALSE, msLevel = 1)  ## MS 1 data.
ondisk <- readMSData2(f, verbose = FALSE)
ondisk1 <- readMSData2(f, msLevel. = 1, verbose = FALSE)
ondisk2 <- readMSData2(f, msLevel. = 2, verbose = FALSE)


############################################################
## trimMz
test_that("OnDiskMSnExp filter/trimMz", {
    inMem <- inmem1
    onDisk <- ondisk1
    inMemMz <- filterMz(inMem, mzlim = c(500, 550))
    onDiskMz <- filterMz(onDisk, mzlim = c(500, 550))
    expect_true(all.equal(inMemMz, onDiskMz))
})

############################################################
## normalize
test_that("Compare OnDiskMSnExp and MSnExp normalize", {
    ## Compare timings and results for normalize.
    inMemN <- normalize(inMem)
    onDiskN <- normalize(onDisk)
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
    expect_warning(inmem1 <- smooth(inmem1))
    expect_warning(inmem2 <- smooth(inmem2, verbose = FALSE))
    expect_warning(expect_true(all.equal(filterMsLevel(ondisk, msLevel. = 1),
                                         inmem1)))
    expect_warning(expect_true(all.equal(filterMsLevel(ondisk, msLevel. = 2),
                                         inmem2)))
})

############################################################
## compareSpectra
test_that("Compare OnDiskMSnExp and MSnExp compareSpectra", {
    csp <- compareSpectra(inmem1)
    csp2 <- compareSpectra(ondisk1)
    rownames(csp) <- colnames(csp) <- NULL
    rownames(csp2) <- colnames(csp2) <- NULL
    expect_identical(csp, csp2)
})

############################################################
## pickPeaks
test_that("Compare OnDiskMSnExp and MSnExp pickPeaks", {
    ## Setting centroided FALSE to avoid warnings...
    centroided(inmem1) <- FALSE
    pp <- pickPeaks(inmem1)
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
    en <- estimateNoise(inmem1)
    centroided(ondisk1) <- FALSE
    en2 <- estimateNoise(ondisk1)
    names(en) <- NULL
    names(en2) <- NULL
    expect_identical(en, en2)
    ## Same with method = SuperSmoother
    centroided(inmem1) <- FALSE
    en <- estimateNoise(inmem1, method = "SuperSmoother")
    centroided(ondisk1) <- FALSE
    en2 <- estimateNoise(ondisk1, method = "SuperSmoother")
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

############################################################
## bpi
test_that("bpi,OnDiskMSnExp", {
    ## Get the "initial" one.
    basepi <- bpi(ondisk)
    expect_identical(unname(basepi), fData(ondisk)$basePeakIntensity)
    expect_identical(names(basepi), featureNames(ondisk))
    ## Calculate:
    basepi <- bpi(ondisk, initial = FALSE)
    expect_identical(names(basepi), featureNames(ondisk))
    sp <- spectra(ondisk)
    basepi_calc <- unlist(lapply(sp, FUN = function(z) max(intensity(z))))
    expect_identical(basepi, basepi_calc)
})


############################################################
## tic
test_that("tic,OnDiskMSnExp", {
    ## Get the "initial" one.
    totic <- tic(ondisk)
    expect_identical(unname(totic), fData(ondisk)$totIonCurrent)
    expect_identical(names(totic), featureNames(ondisk))
    ## Calculate:
    totic <- tic(ondisk, initial = FALSE)
    expect_identical(names(totic), featureNames(ondisk))
    sp <- spectra(ondisk)
    totic_calc <- unlist(lapply(sp, FUN = tic))
    expect_identical(totic, totic_calc)
})
