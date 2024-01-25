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
    pd <- new("AnnotatedDataFrame",
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
                equals(c("fileIdx", "retention.time",
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
    k <- k1 <- featureNames(aa)[1]
    k2 <- featureNames(aa)[2]
    expect_true(all.equal(aa[[k]], assayData(aa)[[k]]))
    sub.aa <- aa[1:2]
    expect_true(all.equal(sub.aa[[k1]], assayData(sub.aa)[[k1]]))
    expect_true(all.equal(sub.aa[[k2]], assayData(sub.aa)[[k2]]))
    expect_equal(fData(sub.aa), fData(aa)[1:2, , drop = FALSE])
    my.prec <- precursorMz(aa)[1]
    my.prec.aa <- extractPrecSpectra(aa, my.prec)
    expect_true(all(precursorMz(my.prec.aa) == my.prec))
    expect_equal(length(my.prec.aa), 2)
    expect_equal(ls(assayData(my.prec.aa)), paste0("F1.S", c(1,3)))
    ## subsetting errors
    expect_error(aa[[1:3]], "subscript out of bounds")
    expect_error(aa[c("F1.S1","F1.S2")],
                 "subsetting works only with numeric or logical")
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
    aa <- extdata_mzXML_in_mem_ms2
    clean.aa <- clean(aa, verbose = FALSE)
    rmpeaks.aa <- removePeaks(aa, verbose = FALSE)
    expect_equal(ls(assayData(clean.aa)), ls(assayData(aa)))
    expect_equal(ls(assayData(rmpeaks.aa)), ls(assayData(aa)))
    int <- c(0, 2, 3, 1, 0, 0, 1)
    sp <- new("Spectrum2",
              intensity = int,
              mz = 1:length(int),
              centroided = FALSE)
    rsp <- removePeaks(sp)
    expect_equal(peaksCount(sp), length(int))
    expect_equal(ionCount(sp), sum(int))
    expect_equal(removePeaks(sp), rsp)
    expect_equal(ionCount(removePeaks(sp, 1)), 6)
    expect_equal(ionCount(removePeaks(sp, 3)), 0)
    expect_equal(ionCount(removePeaks(sp, max(intensity(sp)))), 0)
    expect_equal(peaksCount(sp), peaksCount(rsp))
    expect_equal(peaksCount(clean(rsp)), 5)
    expect_equal(peaksCount(clean(sp)), 7)
    expect_equal(removePeaks(sp, 0), sp)
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
    expect_error(addIdentificationData(new("MSnExp"), identFile),
                 "No feature data found.")
    aa <- extdata_mzXML_in_mem_ms2
    expect_error(addIdentificationData(aa, "foobar.mzid"),
                 "not found")
    fd <- fData(addIdentificationData(aa, identFile, verbose = FALSE))
    expect_equal(fd$spectrum, 1:5)
    expect_equal(fd$sequence,
                 c("VESITARHGEVLQLRPK", "IDGQWVTHQWLKK",
                   NA, NA, "LVILLFR"))
    expect_equal(fd$DatabaseAccess,
                 c("ECA0984", "ECA1028", NA, NA, "ECA0510"))
    expect_equal(fd$idFile, c("dummyiTRAQ.mzid", "dummyiTRAQ.mzid", NA, NA,
                              "dummyiTRAQ.mzid"))
    expect_equal(fd$npsm.prot, c(1, 1, NA, NA, 1))
    expect_equal(fd$npep.prot, c(1, 1, NA, NA, 1))
    expect_equal(fd$nprot, c(1, 1, NA, NA, 1))
    expect_equal(fd$npsm.pep, c(1, 1, NA, NA, 1))
})


test_that("addIdentificationData to OnDiskMSnExp", {
    quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
    rw1 <- extdata_mzXML_in_mem_ms2
    rw2 <- extdata_mzXML_on_disk_ms2
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
    identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
    aa <- extdata_mzXML_in_mem_ms2
    bb <- addIdentificationData(aa, identFile)
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
    f <- msdata::proteomics(full.names = TRUE,
                            pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    i1 <- isolationWindow(f, unique = FALSE)
    i2 <- isolationWindow(tmt_erwinia_on_disk, unique = FALSE)
    i3 <- isolationWindow(tmt_erwinia_in_mem_ms2, unique = FALSE)
    expect_identical(i1, i2)
    expect_identical(i1, i3)
})

test_that("spectrapply,MSnExp", {
    library(msdata)
    inMem <- microtofq_in_mem_ms1

    sps <- spectra(inMem)
    sps_2 <- spectrapply(inMem)
    expect_identical(sps, sps_2)
    ## apply a function.
    dfs <- spectrapply(inMem, FUN = as, Class = "data.frame")
    dfs_2 <- lapply(sps, FUN = as, Class = "data.frame")
    expect_identical(dfs, dfs_2)
})

test_that("splitByFile,MSnExp", {
    library(msdata)
    inMem <- microtofq_in_mem_ms1
    expect_error(splitByFile(inMem, f = factor(1:3)))
    spl <- splitByFile(inMem, f = factor(c("b", "a")))
    expect_equal(spectra(spl[[1]]), spectra(filterFile(inMem, 2)))
    expect_equal(pData(spl[[1]]), pData(filterFile(inMem, 2)))
    expect_equal(spectra(spl[[2]]), spectra(filterFile(inMem, 1)))
    expect_equal(pData(spl[[2]]), pData(filterFile(inMem, 1)))
})

test_that("$ operator on MSnExp works", {
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    expect_equal(pData(msx)$sampleNames, msx$sampleNames)
    ## replace.
    msx$sampleNames <- "b"
    expect_equal("b", msx$sampleNames)
    expect_equal(msx$bla, NULL)
    ## Add a new column
    msx$newCol <- 5
    expect_equal(msx$newCol, 5)
})

test_that("pData<- on MSnExp works", {
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx

    newDf <- data.frame(sampleName = "b", otherCol = 3)
    ## replace.
    pData(msx) <- newDf
    expect_equal(pData(msx), newDf)
    expect_error(pData(msx) <- 13)
})

test_that("phenoData<- on MSnExp works", {
    im <- microtofq_in_mem_ms1
    old_pd <- phenoData(im)

    expect_error(phenoData(im) <- 4)
    ## phenoData(im) <- data.frame(a = 4)
    pd_2 <- old_pd
    pData(pd_2) <- cbind(pData(old_pd), add_col = 4)
    ## Assign AnnotatedDataFrame
    phenoData(im) <- AnnotatedDataFrame(pData(pd_2))
    expect_true(is(im@phenoData, "AnnotatedDataFrame"))
    expect_equal(pData(im), pData(pd_2))
    ## Assign AnnotatedDataFrame
    phenoData(im) <- pd_2
    expect_true(is(im@phenoData, "AnnotatedDataFrame"))
    expect_equal(phenoData(im), pd_2)
})

test_that("chromatogram,MSnExp works", {
    inMem <- microtofq_in_mem_ms1
    ## Reduce here the tests. Most of the tests are performed in
    ## chromatogram,OnDiskMSnExp and both methods use the same low level
    ## function.

    ## Multiple mz ranges.
    mzr <- matrix(c(100, 120, 200, 220, 300, 320), nrow = 3, byrow = TRUE)
    rtr <- matrix(c(50, 300), nrow = 1)
    res <- chromatogram(inMem, mz = mzr, rt = rtr)
    ## Check that the values for all ranges is within the specified ranges
    for (i in 1:nrow(mzr)) {
        expect_true(all(mz(res[i, 1]) >= mzr[i, 1] &
                        mz(res[i, 1]) <= mzr[i, 2]))
        expect_true(all(mz(res[i, 2]) >= mzr[i, 1] &
                        mz(res[i, 2]) <= mzr[i, 2]))
        expect_true(all(rtime(res[i, 1]) >= rtr[1, 1] &
                        rtime(res[i, 1]) <= rtr[1, 2]))
        expect_true(all(rtime(res[i, 2]) >= rtr[1, 1] &
                        rtime(res[i, 2]) <= rtr[1, 2]))
    }
    ## Check that values are correct.
    flt <- filterMz(filterRt(inMem, rt = rtr[1, ]), mz = mzr[2, ])
    ints <- split(unlist(lapply(spectra(flt), function(z) sum(intensity(z)))),
                  fromFile(flt))
    expect_equal(ints[[1]], intensity(res[2, 1]))
    expect_equal(ints[[2]], intensity(res[2, 2]))
    expect_equal(split(rtime(flt), fromFile(flt))[[1]], rtime(res[2, 1]))
    expect_equal(split(rtime(flt), fromFile(flt))[[2]], rtime(res[2, 2]))
    ## fData
    expect_true(nrow(fData(res)) == nrow(res))
    expect_true(all(colnames(fData(res)) == c("mzmin", "mzmax",
                                              "rtmin", "rtmax", "polarity")))
    expect_true(all(fData(res)$rtmin == 50))
    expect_true(all(fData(res)$rtmax == 300))
    expect_equal(fData(res)$mzmin, c(100, 200, 300))
    expect_equal(fData(res)$mzmax, c(120, 220, 320))
    expect_equal(fData(res)$polarity, c(1, 1, 1))

    ## Now with ranges for which we don't have values in one or the other.
    rtr <- matrix(c(280, 300, 20, 40), nrow = 2,
                  byrow = TRUE)  ## Only present in first, or 2nd file
    res <- chromatogram(inMem, rt = rtr)
    ## Check fromFile
    for (i in 1:ncol(res))
        expect_true(all(sapply(res[, i], fromFile) == i))
    expect_equal(length(res[2, 1]), 0)
    expect_equal(length(res[1, 2]), 0)
    ## Check rtime
    expect_true(all(rtime(res[1, 1]) >= rtr[1, 1] &
                    rtime(res[1, 1]) <= rtr[1, 2]))
    expect_true(all(rtime(res[2, 2]) >= rtr[2, 1] &
                    rtime(res[2, 2]) <= rtr[2, 2]))
    ## Check intensity
    flt <- filterRt(inMem, rt = rtr[1, ])
    spctr <- split(spectra(flt), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[1, 1]))
    flt <- filterRt(inMem, rt = rtr[2, ])
    spctr <- split(spectra(flt), fromFile(flt))
    ints <- unlist(lapply(spctr[[1]], function(z) sum(intensity(z))))
    expect_equal(ints, intensity(res[2, 2]))

    ## Check that phenoType is correctly passed.
    pd <- data.frame(name = c("first", "second"), idx = 1:2)
    pData(inMem) <- pd
    chrs <- chromatogram(inMem)
    ## rownames(pd) <- colnames(chrs)
    expect_equal(pData(chrs), pd)

    chrs_2 <- chromatogram(inMem, msLevel = 1:4)
    expect_equal(chrs, chrs_2)
})

test_that("setAs,MSnExp,data.frame works", {
    od <- microtofq_on_disk
    im <- microtofq_in_mem_ms1

    res <- filterMz(im, mz = c(200, 300))
    df <- as(res, "data.frame")
    expect_equal(colnames(df), c("file", "rt", "mz", "i"))
    expect_equal(unlist(mz(res), use.names = FALSE), df$mz)
    expect_equal(unlist(intensity(res), use.names = FALSE), df$i)
})

test_that("pickPeaks,MSnExp works with msLevel", {
    res <- pickPeaks(tmt_im_ms1_sub, msLevel = 2)
    expect_identical(peaksCount(res), peaksCount(tmt_im_ms1_sub))
    res <- pickPeaks(tmt_im_ms1_sub, msLevel = 1:3)
    expect_true(all(peaksCount(res) < peaksCount(tmt_im_ms1_sub)))
})

test_that("smooth,MSnExp works with msLevel", {
    res <- smooth(tmt_im_ms1_sub, msLevel = 2)
    expect_identical(intensity(res), intensity(tmt_im_ms1_sub))
    res <- smooth(tmt_im_ms1_sub, msLevel = 1:3)
    expect_false(all(unlist(intensity(res)) ==
                     unlist(intensity(tmt_im_ms1_sub))))
})

test_that("pickPeaks,MSnExp works with refineMz", {
    ## Reduce the TMT erwinia data set to speed up processing on the full
    ## data.
    tmt <- tmt_im_ms1_sub
    centroided(tmt) <- FALSE
    ## Get one spectrum against which we will compare
    spctr <- tmt[[1]]
    ## kNeighbors
    tmt_pk <- pickPeaks(tmt, refineMz = "kNeighbors", k = 1)
    spctr_pk <- pickPeaks(spctr, refineMz = "kNeighbors", k = 1)
    spctr_tmt_pk <- tmt_pk[[1]]
    expect_equal(spctr_pk, spctr_tmt_pk)
    ## descendPeak
    spctr_pk <- pickPeaks(spctr, refineMz = "descendPeak",
                          signalPercentage = 75)
    tmt_pk <- pickPeaks(tmt, refineMz = "descendPeak",
                        signalPercentage = 75)
    expect_equal(spctr_pk, tmt_pk[[1]])
    ## Check if we can call method and refineMz and pass arguments to both
    spctr_pk <- pickPeaks(spctr, refineMz = "kNeighbors", k = 1,
                          method = "SuperSmoother", span = 0.9)
    tmt_pk <- pickPeaks(tmt, refineMz = "kNeighbors",
                        k = 1, method = "SuperSmoother", span = 0.9)
    expect_equal(spctr_pk, tmt_pk[[1]])
    ## Check errors
    expect_error(pickPeaks(tmt_erwinia_in_mem_ms1, refineMz = "some_method"))
})

test_that("estimateMzResolution,MSnExp works", {
    res <- estimateMzResolution(tmt_erwinia_in_mem_ms2)
    expect_equal(unname(res[[15]]),
                 estimateMzResolution(tmt_erwinia_in_mem_ms2[[15]]))
})

test_that("estimateMzScattering works", {
    expect_error(estimateMzScattering(4))

    res <- estimateMzScattering(tmt_erwinia_in_mem_ms1)
    mzr <- estimateMzResolution(tmt_erwinia_in_mem_ms1)
    idx <- which.max(spectrapply(tmt_erwinia_in_mem_ms1, peaksCount))
    ## m/z scattering should be smaller than m/z resolution
    expect_true(res[[idx]] < mzr[[idx]])

    ## .estimate_mz_scattering_list.
    res_2 <- .estimate_mz_scattering_list(spectra(tmt_erwinia_in_mem_ms1),
                                          timeDomain = FALSE)
    expect_equal(unname(res), res_2)
})

test_that("combineSpectraMovingWindow works", {
    ## Check errors
    expect_error(combineSpectraMovingWindow("3"))

    od <- filterFile(sciex, 1)
    ## Focus on the one with most peaks
    idx <- which.max(peaksCount(od))
    spctrl <- spectra(od[(idx -1):(idx + 1)])

    od <- od[(idx - 3):(idx + 3)]
    spctr <- meanMzInts(spctrl, timeDomain = TRUE, unionPeaks = FALSE,
                        main = 2L)
    ## Should be different from raw ones
    expect_true(is.character(all.equal(mz(spctr), mz(od[[4]]))))
    expect_true(is.character(all.equal(intensity(spctr), intensity(od[[4]]))))

    ## Use pre-calculated mzd:
    mzd <- estimateMzScattering(od)
    ## If mzd is estimated on mz and combination on sqrt(mz) it will fail.
    od_comb <- combineSpectraMovingWindow(
        od, mzd = mzd[[4]], timeDomain = TRUE)
    ## All on m/z scale
    od_comb <- combineSpectraMovingWindow(od, mzd = mzd[[4]],
                                          timeDomain = FALSE)
    spctr_comb <- od_comb[[4]]
    expect_equal(mz(spctr_comb), mz(spctr))
    expect_equal(intensity(spctr_comb), intensity(spctr))

    ## Estimate on the sqrt(mz)
    mzd <- estimateMzScattering(od, timeDomain = TRUE)
    od_comb <- combineSpectraMovingWindow(od, mzd = mzd[[4]], timeDomain = TRUE)
    spctr_comb <- od_comb[[4]]
    expect_equal(mz(spctr_comb), mz(spctr))
    expect_equal(intensity(spctr_comb), intensity(spctr))
    ## Should be different from the original ones.
    spctr_raw <- od[[4]]
    expect_true(is.character(all.equal(mz(spctr), mz(spctr_raw))))
    expect_true(is.character(all.equal(intensity(spctr), intensity(spctr_raw))))

    ## Shouldn't make a difference if we're using timeDomain = TRUE or FALSE.
    od_comb <- combineSpectraMovingWindow(od)
    spctr_comb <- od_comb[[4]]
    expect_equal(mz(spctr_comb), mz(spctr))
    expect_equal(intensity(spctr_comb), intensity(spctr))

    od_comb_w <- combineSpectraMovingWindow(od, weighted = TRUE)
    spctr_w <- od_comb_w[[4]]
    expect_true(sum(mz(spctr_comb) != mz(spctr_w)) > length(mz(spctr_w))/2)
    expect_equal(intensity(spctr_comb), intensity(spctr_w))

    expect_equal(length(od), length(od_comb))
    expect_equal(peaksCount(od), peaksCount(od_comb))
})


test_that("plotXIC_MSnExp works", {
    im <- microtofq_in_mem_ms1
    expect_warning(plotXIC_MSnExp(filterMz(im, c(600, 680))))
    plotXIC_MSnExp(filterMz(im, c(610, 615)), pch = 23)
    ## filter to get only one
    plotXIC_MSnExp(filterMz(filterRt(im, c(270, 280)), c(610, 615)), cex = 2)

    expect_error(plotXIC_MSnExp(tmt_erwinia_in_mem_ms2))
})

test_that("as,MSnExp,MSpectra works", {
    res <- as(tmt_erwinia_in_mem_ms1, "MSpectra")
    expect_equal(res@listData, spectra(tmt_erwinia_in_mem_ms1))
    expect_true(ncol(mcols(res)) == 0)

    res <- as(sciex, "MSpectra")
    expect_equal(length(res), length(sciex))
    expect_equal(msLevel(res), msLevel(sciex))
    expect_equal(intensity(res), intensity(sciex))
    expect_true(ncol(mcols(res)) > 0)
})

test_that("isolationWindowLowerMz, isolationWindowUpperMz work", {
    expect_error(isolationWindowLowerMz(tmt_im_ms2_sub), "not available")
    expect_error(isolationWindowUpperMz(tmt_im_ms2_sub), "not available")
})

test_that("combineSpectra,MSnExp works", {
    res <- combineSpectra(tmt_im_ms2_sub)
    expect_true(is(res, "MSnExp"))
    expect_true(length(res) == 1)
    expect_equal(res[[1]],
                 combineSpectra(as(tmt_im_ms2_sub, "MSpectra"))@listData[[1]])
    res2 <- combineSpectra(tmt_im_ms2_sub, mzd = 0.1)
    expect_true(peaksCount(res2) < peaksCount(res))
    res3 <- combineSpectra(tmt_im_ms2_sub, mzd = 0.01, minProp = 0.1,
                           method = consensusSpectrum)
    expect_true(peaksCount(res3) < peaksCount(res2))

    expect_error(combineSpectra(tmt_im_ms2_sub, fcol = "other"))
})
