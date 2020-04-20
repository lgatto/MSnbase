context("readWriteMgfData")

test_that("extractMgfSpectrum2Info", {
    mgf <- c("TITLE=foobar File=\"foobar.raw\", Native ID:\"controllerType=0 controllerNumber=1, scan=100\"",
             "RTINSECONDS=600",
             "PEPMASS=100 50000",
             "CHARGE=3+",
             "10 100",
             "11 200",
             "12 300",
             "13 400")

    s <- new("Spectrum2",
             rt = 600,
             scanIndex = 0L,
             precursorMz = 100,
             precursorIntensity = 50000,
             precursorCharge = 3L,
             mz = c(10, 11, 12, 13),
             intensity = c(100, 200, 300, 400),
             fromFile = 1L,
             acquisitionNum = NA_integer_,
             collisionEnergy = NA_real_,
             precScanNum = NA_integer_,
             centroided = TRUE)

    fdata <- c(TITLE="foobar File=\"foobar.raw\", Native ID:\"controllerType=0 controllerNumber=1, scan=100\"",
               RTINSECONDS="600",
               PEPMASS="100 50000",
               CHARGE="3+")
    result <- list(spectrum=s, fdata=fdata)
    expect_equal(MSnbase:::extractMgfSpectrum2Info(mgf, centroided = TRUE),
                 result)
})

test_that("writeMgfContent works", {
    s <- new("Spectrum2",
             rt = 600,
             scanIndex = 0L,
             precursorMz = 100,
             precursorIntensity = 50000,
             precursorCharge = 3L,
             mz = c(10, 11, 12, 13),
             intensity = c(100, 200, 300, 400),
             fromFile = 1L,
             acquisitionNum = NA_integer_,
             collisionEnergy = NA_real_,
             precScanNum = NA_integer_,
             centroided = TRUE)
    tmpf <- tempfile()
    MSnbase:::writeMgfContent(s, con = tmpf)
    lns <- readLines(tmpf)
    expect_equal(lns[grep("^RT", lns)], "RTINSECONDS=600")
    expect_equal(lns[8], "10 100")
    res <- readMgfData(tmpf)
    expect_equal(mz(res[[1]]), mz(s))
    
    file.remove(tmpf)
    MSnbase:::writeMgfContent(s, con = tmpf, addFields = c(ID = 12))
    lns <- readLines(tmpf)
    expect_equal(lns[8], "ID=12")

    file.remove(tmpf)
    MSnbase:::writeMgfContent(s, con = tmpf, addFields = c(ID = 12,
                                                           feature = "B"))
    lns <- readLines(tmpf)
    expect_equal(lns[8], "ID=12")
    expect_equal(lns[9], "FEATURE=B")

    res <- readMgfData(tmpf)
    expect_equal(mz(res[[1]]), mz(s))
    expect_equal(intensity(res[[1]]), intensity(s))
})

test_that("writeMgfDataFile works", {

    s1 <- new("Spectrum2",
              rt = 600,
              scanIndex = 0L,
              precursorMz = 100,
              precursorIntensity = 50000,
              precursorCharge = 3L,
              mz = c(10, 11, 12, 13),
              intensity = c(100, 200, 300, 400),
              fromFile = 1L,
              acquisitionNum = NA_integer_,
              collisionEnergy = NA_real_,
              precScanNum = NA_integer_,
              centroided = TRUE)
    s2 <- new("Spectrum2",
              rt = 601, mz = c(12, 13, 14, 15, 16),
              intensity = c(200, 3000, 4000, 300, 399),
              fromFile = 2L)

    tmpf <- tempfile()
    MSnbase:::writeMgfDataFile(list(s1, s2), tmpf)
    lns <- readLines(tmpf)
    
    ## With addFields
    file.remove(tmpf)
    af <- data.frame(pkid = c("a", "b"), value = 1:2)

    expect_error(MSnbase:::writeMgfDataFile(list(s1, s2), tmpf, addFields = 3))
    expect_error(MSnbase:::writeMgfDataFile(list(s1, s2), tmpf,
                                            addFields = af[1, ]))
    MSnbase:::writeMgfDataFile(list(s1, s2), tmpf, addFields = af)
    lns <- readLines(tmpf)
    expect_equal(lns[8], "PKID=a")
    expect_equal(lns[9], "VALUE=1")
    expect_equal(lns[21], "PKID=b")
    expect_equal(lns[22], "VALUE=2")
})


test_that("writeMgfContent works with file and con", {
    s <- new("Spectrum2",
             rt = 600,
             scanIndex = 0L,
             precursorMz = 100,
             precursorIntensity = 50000,
             precursorCharge = 3L,
             mz = c(10, 11, 12, 13),
             intensity = c(100, 200, 300, 400),
             fromFile = 1L,
             acquisitionNum = NA_integer_,
             collisionEnergy = NA_real_,
             precScanNum = NA_integer_,
             centroided = TRUE)
    
    tmpf1 <- tempfile()
    tmpf2 <- tempfile()
    tmpcon <- file(description = tmpf2,
                   open = "w")

    MSnbase:::writeMgfContent(s, con = tmpf1)
    MSnbase:::writeMgfContent(s, con = tmpcon)

    x1 <- readMgfData(tmpf1)
    x2 <- readMgfData(tmpf2)
    
    expect_identical(fileNames(x1), tmpf1)
    expect_identical(fileNames(x2), tmpf2)
    
    ## pData and processingData are expected to be different
    x1$sampleNames <- x2$sampleNames <- ""
    x1@processingData@files <- x2@processingData@files <- ""
    x1@processingData@processing <- x2@processingData@processing <- "" 

    expect_equal(x1, x2)    
})
