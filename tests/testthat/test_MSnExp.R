context("MSnExp class")

test_that("MSnExp validity", {
    expect_true(validObject(new("MSnExp")))
    data(itraqdata, package = "MSnbase")
    expect_true(validObject(itraqdata))
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    expect_true(validObject(msx))
})

test_that("readMSData", {
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    file <- dir(system.file(package = "MSnbase", dir = "extdata"),
                full.name = TRUE,pattern = "mzXML$")
    aa <- readMSData(file, verbose = FALSE, centroided. = FALSE)
    expect_identical(as.list(assayData(aa)), as.list(assayData(msx)))
    ## ## removing below due to spurious error on windows
    ## ## processingData will be different
    ## aa@processingData <- processingData(msx)
    ## ## overwrite R/Bioc versions
    ## msx@.__classVersion__ <- aa@.__classVersion__
    ## ## msx has ident data to be remove for comparison
    ## fData(msx) <- fData(msx)[, 1, drop = FALSE]
    ## expect_true(all.equal(aa, msx))
})

test_that("readMSData with pdata", {
    file <- dir(system.file(package = "MSnbase", dir = "extdata"),
                full.name = TRUE, pattern = "mzXML$")
    pd <- new("NAnnotatedDataFrame",
              data = data.frame(pvarA = "A", pvarB = "B"))
    aa <- readMSData(file, pdata = pd, verbose = FALSE)
    expect_true(validObject(aa))
    expect_true(validObject(pData(aa)))
    expect_true(all.equal(dim(pd), dim(phenoData(aa))))
})

test_that("readMSData and dummy MSnExp msLevel 2 instance", {
    file <- dir(system.file(package = "MSnbase", dir = "extdata"),
                full.name = TRUE, pattern = "mzXML$")
    aa <- readMSData(file, verbose = FALSE, centroided. = FALSE)
    expect_true(class(aa) == "MSnExp")
    ## centroided get and set
    expect_false(any(centroided(aa)))
    val <- rep(TRUE, length(aa))
    centroided(aa) <- val
    expect_true(validObject(aa))
    expect_true(all(centroided(aa)))
    val[sample(length(aa), 2)] <- FALSE
    centroided(aa) <- val
    expect_true(sum(centroided(aa)) == length(aa)-2)
    centroided(aa) <- rep(FALSE, length(aa))
    expect_false(any(centroided(aa)))
    ## checking slots and methods
    expect_equal(length(aa), 5)
    expect_that(nrow(header(aa)), equals(length(aa)))
    expect_that(names(header(aa)),
                equals(c("file", "retention.time",
                         "precursor.mz", "precursor.intensity",
                         "charge", "peaks.count","tic","ionCount",
                         "ms.level", "acquisition.number",
                         "collision.energy")))
    ## MS levels
    expect_equal(length(msLevel(aa)), 5)
    expect_equal(unique(msLevel(aa)), 2)
    expect_equal(length(precScanNum(aa)), 5)
    expect_equal(length(unique(precScanNum(aa))), 1)
    ## Precursor MZ
    expect_equal(length(precursorMz(aa)), 5)
    expect_that(precursorMz(aa)[1], is_a("numeric"))
    expect_equal(length(unique(precursorMz(aa))), 4)
    expect_equal(range(precursorMz(aa)),
                 c(437.80401611, 716.34051514))
    expect_equal(as.numeric(sort(precursorMz(aa))[1]), 437.80401611) ## [*]
    ## Retention time
    expect_equal(length(rtime(aa)), 5)
    expect_that(rtime(aa)[1], is_a("numeric"))
    expect_equal(range(rtime(aa)), c(1501.35, 1502.31))
    expect_equal(as.numeric(sort(rtime(aa))[1]), 1501.35) ## [*]
    ## [*] using as.numeric because rtime and precursorMz return named numerics
    ## Meta data
    expect_equal(dim(fData(aa)), c(5, 1))
    expect_equal(dim(pData(aa)), c(1, 1))
    ## subsetting
    expect_true(all.equal(aa[["X4.1"]], assayData(aa)[["X4.1"]]))
    sub.aa <- aa[1:2]
    expect_true(all.equal(sub.aa[["X1.1"]], assayData(sub.aa)[["X1.1"]]))
    expect_true(all.equal(sub.aa[["X2.1"]], assayData(sub.aa)[["X2.1"]]))
    expect_equal(fData(sub.aa), fData(aa)[1:2, , drop = FALSE])
    my.prec <- precursorMz(aa)[1]
    my.prec.aa <- extractPrecSpectra(aa, my.prec)
    expect_true(all(precursorMz(my.prec.aa) == my.prec))
    expect_equal(length(my.prec.aa), 2)
    expect_equal(ls(assayData(my.prec.aa)), paste0("X", c(1,3), ".1"))
    ## subsetting errors
    expect_error(aa[[1:3]], "subscript out of bounds")
    expect_error(aa[c("X1.1","X2.1")], "subsetting works only with numeric or logical")
    expect_error(aa[["AA"]], "object 'AA' not found")
    expect_error(aa[1:10], "subscript out of bounds")
    ## testing that accessors return always attributes in same order
    precMzNames <- names(precursorMz(aa))
    ionCountNames <- names(ionCount(aa))
    expect_true(all.equal(precMzNames, ionCountNames))
    precChNames <- names(precursorCharge(aa))
    expect_that(precMzNames, equals(precChNames))
    aqnNames <- names(acquisitionNum(aa))
    expect_that(precMzNames, equals(aqnNames))
    rtNames <- names(rtime(aa))
    expect_that(precMzNames, equals(rtNames))
    pkCntNames <- names(peaksCount(aa))
    expect_that(precMzNames, equals(pkCntNames))
    mslNames <- names(msLevel(aa))
    expect_that(precMzNames, equals(mslNames))
    coleNames <- names(collisionEnergy(aa))
    expect_that(precMzNames, equals(coleNames))
    intNames <- names(intensity(aa))
    expect_that(precMzNames, equals(intNames))
    mzNames <- names(mz(aa))
    expect_that(precMzNames, equals(mzNames))
    ffNames <- names(fromFile(aa))
    expect_that(precMzNames, equals(ffNames))
})

context("MSnExp processing")

## ! Issues with edited dummy data for MS1 uploading, although
## ! things work fine for original data set. Commented these
## ! tests for the moment
## test_that("readMSLData and dummy MSnExp msLevel 1 instance", {
##   file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
##   aa <- readMSData(file,msLevel=1,verbose=FALSE)
##   expect_that(class(aa)=="MSnExp",is_true())
##   expect_equal(length(aa),equals(5))
##   ## MS levels
##   expect_that(length(msLevel(aa)),equals(5))
##   expect_that(unique(msLevel(aa)),equals(1))
##   ## Retention time
##   expect_that(length(rtime(aa)),equals(5))
##   expect_that(rtime(aa)[1],is_a("numeric"))
##   expect_that(range(rtime(aa)),
##               equals(c(1982.08,3015.47)))
##   expect_that(as.numeric(polarity(aa)),equals(rep(-1,length(aa)))) ## [*]
##   expect_that(as.numeric(rtime(aa)[1]),equals(1982.08)) ## [*]
##   ## [*] using as.numeric because rtime and precursorMz return named numerics
## })

context("MSnExp data")

test_that("spectra order and integrity", {
    file <- dir(system.file(package = "MSnbase", dir = "extdata"),
                full.name = TRUE,
                pattern = "mzXML$")
    aa <- readMSData(file, verbose = FALSE, centroided. = FALSE)
    clean.aa <- clean(aa, verbose = FALSE)
    rmpeaks.aa <- removePeaks(aa, verbose = FALSE)
    expect_that(ls(assayData(clean.aa)), equals(ls(assayData(aa))))
    expect_that(ls(assayData(rmpeaks.aa)), equals(ls(assayData(aa))))
    int <- c(0, 2, 3, 1, 0, 0, 1)
    sp <- new("Spectrum2",
              intensity = int,
              mz = 1:length(int),
              centroided = FALSE)
    rsp <- removePeaks(sp)
    expect_that(peaksCount(sp), equals(length(int)))
    expect_that(ionCount(sp), equals(sum(int)))
    expect_that(all.equal(removePeaks(sp),rsp), is_true())
    expect_that(ionCount(removePeaks(sp,1)), equals(6))
    expect_that(ionCount(removePeaks(sp,3)), equals(0))
    expect_that(ionCount(removePeaks(sp,max(intensity(sp)))), equals(0))
    expect_that(peaksCount(sp), equals(peaksCount(rsp)))
    expect_that(peaksCount(clean(rsp)), equals(6))
    expect_that(peaksCount(clean(sp)), equals(7))
    expect_that(all.equal(removePeaks(sp,0), sp), is_true())
})

test_that("MSnExp normalisation", {
    data(itraqdata, package = "MSnbase")
    aa <- itraqdata[1:3]
    bb <- normalise(aa, "max")
    expect_true(all(sapply(intensity(bb), max) == 1))
    expect_true(all.equal(unlist(sapply(intensity(aa), order)),
                          unlist(sapply(intensity(bb), order))))
})

context("MSnExp identification data")

test_that("addIdentificationData", {
    quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")

    expect_error(addIdentificationData(new("MSnExp"),
                                       identFile, verbose = FALSE),
                 "No feature data found.")

    aa <- readMSData(quantFile, verbose = FALSE)

    expect_error(addIdentificationData(aa, "foobar.mzid",
                                       verbose = FALSE),
                 "does not exist")

    fd <- fData(addIdentificationData(aa, identFile, verbose = FALSE))

    expect_equal(fd$spectrum, 1:5)
    expect_equal(fd$pepseq,
                 c("VESITARHGEVLQLRPK", "IDGQWVTHQWLKK",
                   NA, NA, "LVILLFR"))
    expect_equal(fd$accession,
                 c("ECA0984;ECA3829", "ECA1028",
                   NA, NA, "ECA0510"))
    expect_equal(fd$idFile, c("dummyiTRAQ.mzid", "dummyiTRAQ.mzid", NA, NA,
                              "dummyiTRAQ.mzid"))
    expect_equal(fd$npsm.prot, c(1, 1, NA, NA, 1))
    expect_equal(fd$npep.prot, c(1, 1, NA, NA, 1))
    expect_equal(fd$nprot, c(2, 1, NA, NA, 1))
    expect_equal(fd$npsm.pep, c(1, 1, NA, NA, 1))
})


test_that("addIdentificationData to OnDiskMSnExp", {
    quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
    rw1 <- readMSData(quantFile)
    rw2 <- readMSData2(quantFile)
    expect_true(all.equal(rw1, rw2))
    rw1 <- addIdentificationData(rw1, identFile)
    rw2 <- addIdentificationData(rw2, identFile)
    expect_true(all.equal(rw1, rw2))
    k <- intersect(fvarLabels(rw1), fvarLabels(rw2))
    expect_identical(fData(rw1)[ k], fData(rw2)[ k])
})

## test_that("addIdentificationData from MSGF+ and X!TANDEM", {
##     rawFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##                    full.name = TRUE, pattern = "mzXML$")
##     msgfFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
##     xtFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##                   full.name = TRUE, pattern = "dummyiTRAQxt.mzid")
##     x <- addIdentificationData(aa, msgfFile)
##     y <- addIdentificationData(aa, xtFile)
## })


test_that("idSummary", {
    quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")

    aa <- readMSData(quantFile, verbose = FALSE)
    bb <- addIdentificationData(aa, identFile, verbose = FALSE)

    expect_error(idSummary(aa), "No quantification/identification data found")
    expect_equal(idSummary(bb),
                 data.frame(spectrumFile="dummyiTRAQ.mzXML",
                            idFile="dummyiTRAQ.mzid", coverage=0.6,
                            stringsAsFactors=FALSE))
})

test_that("MSnExp sample names", {
    data(itraqdata, package = "MSnbase")
    expect_identical(sampleNames(itraqdata), "1")
    sampleNames(itraqdata) <- "A"
    expect_identical(sampleNames(itraqdata), "A")
})

test_that("MSnExp mulitplexed sample names", {
    data(itraqdata, package = "MSnbase")
    ## this is an iTRAQ4-plex sample; let's update phenoData
    ## accordingly
    pd2 <- new("NAnnotatedDataFrame",
               data = data.frame(sampleNumbers = 1:4,
                                 row.names = paste0("iTRAQ", 1:4)),
               multiplex = 4,
               multiLabels = paste0("iTRAQ", 1:4))
    itraqdata@phenoData <- pd2
    expect_true(validObject(itraqdata))
    sampleNames(itraqdata) <- LETTERS[1:4]
    expect_identical(sampleNames(itraqdata), LETTERS[1:4])
})

test_that("feautre names are correct", {
    data(itraqdata, package = "MSnbase")
    fn0 <- c("X1" , "X10", "X11", "X12", "X13", "X14", "X15", "X16",
             "X17", "X18", "X19", "X2", "X20", "X21", "X22", "X23",
             "X24", "X25", "X26", "X27", "X28", "X29", "X3" , "X30",
             "X31", "X32", "X33", "X34", "X35", "X36", "X37", "X38",
             "X39", "X4" , "X40", "X41", "X42", "X43", "X44", "X45",
             "X46", "X47", "X48", "X49", "X5" , "X50", "X51", "X52",
             "X53", "X54", "X55", "X6" , "X7" , "X8" , "X9")
    expect_identical(featureNames(itraqdata), fn0)
    ## these below are redundant with the validity method, but keeping
    ## here as this validity rule will change in the future, and want
    ## to have a regression test.
    expect_identical(ls(assayData(itraqdata)), fn0)
    expect_identical(featureNames(featureData(itraqdata)), fn0)
})

test_that("Noise estimation MSnExp", {
    expect_identical(estimateNoise(new("MSnExp")), list())
    e <- new.env()
    e$s1 <- new("Spectrum2", mz = 1:5, intensity = c(1:3, 2:1),
                fromFile = 1L, centroided = FALSE)
    e$s2 <- new("Spectrum2", mz = 3, intensity = 3, centroided = TRUE,
                fromFile = 1L)
    fd <- data.frame(x = 1:2,
                     file = rep("1", 2),
                     row.names = c("s1", "s2"))
    pd <- new("MSnProcess", files = "1")
    msx <- new("MSnExp", assayData = e,
               featureData = new("AnnotatedDataFrame", data = fd),
               processingData = pd)
    expect_warning(ns <- estimateNoise(msx))
    expect_true(inherits(ns, "list"))
    expect_identical(length(ns), 2L)
    expect_identical(ns[[1]], estimateNoise(e$s1))
    expect_identical(ns[[2]], suppressWarnings(estimateNoise(e$s2)))
})

test_that("isolation window", {
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    i1 <- isolationWindow(f, unique = FALSE)
    i2 <- isolationWindow(readMSData2(f), unique = FALSE)
    i3 <- isolationWindow(readMSData(f), unique = FALSE)
    expect_identical(i1, i2)
    expect_identical(i1, i3)
})

test_that("spectrapply,MSnExp", {
    library(msdata)
    mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
             system.file("microtofq/MM8.mzML", package = "msdata"))
    inMem <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE)

    sps <- spectra(inMem)
    sps_2 <- spectrapply(inMem)
    expect_identical(sps, sps_2)
    ## apply a function.
    dfs <- spectrapply(inMem, FUN = as, Class = "data.frame")
    dfs_2 <- lapply(sps, FUN = as, Class = "data.frame")
    expect_identical(dfs, dfs_2)
})
