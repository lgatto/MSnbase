library(msdata)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))
inMem <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE)
onDisk <- readMSData2(files = mzf, msLevel. = 1, centroided. = TRUE)

test_that("Subsetting between OnDiskMSnExp and MSnExp - [[", {
    ## Extract individual spectra.
    sp1 <- inMem[[77]]
    sp2 <- onDisk[[77]]
    expect_identical(sp1, sp2)
})

test_that("Subsetting between OnDiskMSnExp and MSnExp - [[ by name", {
    ## by name.
    theN <- featureNames(inMem)[100]
    sp1 <- inMem[[theN]]
    sp2 <- onDisk[[theN]]
    ## ?
    expect_identical(sp1, sp2)
    theN <- 100
    sp1 <- inMem[[theN]]
    sp2 <- onDisk[[theN]]
    expect_identical(sp1, sp2)
})

test_that("Extract multiple with [[ error", {
    expect_error(onDisk[[c(2, 4, 6)]])
    expect_error(inMem[[c(2, 4, 6)]])
})


test_that("Subsetting between OnDiskMSnExp and MSnExp - [", {
    sub1 <- inMem[1:20, ]
    sub2 <- onDisk[1:20, ]
    expect_true(all.equal(sub1, sub2))
})

test_that("Subsetting OnDiskMSnExp and MSnExp - [ and processingData", {
    expect_true(all.equal(inMem[1], onDisk[1]))
    ## That forces sub-setting of processingData etc
    sp1 <- inMem[1]
    sp2 <- onDisk[1]
    expect_identical(experimentData(sp1), experimentData(sp2))
    expect_true(all.equal(sp1, sp2))
    expect_identical(fileNames(sp1), fileNames(sp2))
    expect_identical(fromFile(sp1), fromFile(sp2))
    ## from second file only:
    sp1 <- inMem[c(2, 4, 6)]
    sp2 <- onDisk[c(2, 4, 6)]
    expect_identical(experimentData(sp1), experimentData(sp2))
    expect_true(all.equal(sp1, sp2))
    expect_identical(fileNames(sp1), fileNames(sp2))
    expect_identical(fromFile(sp1), fromFile(sp2))
})

test_that("Subsetting between OnDiskMSnExp and MSnExp - [ and phenoData", {
    ## Some tests evaluating the correct sub-setting of phenoData etc.
    ## Extract spectra from the first file
    subs <- onDisk[c(1, 3, 5)]
    expect_identical(fileNames(subs), fileNames(onDisk)[1])
    expect_true(all(fromFile(subs) == 1))
    expect_identical(pData(subs), droplevels(pData(onDisk)[1, , drop = FALSE]))
    expect_identical(experimentData(subs)@instrumentManufacturer,
                     experimentData(onDisk)@instrumentManufacturer[1])
})

test_that("Subsetting between OnDiskMSnExp and MSnExp - second file", {
    ## Extract spectra from the second file
    subs <- onDisk[c(2, 4, 6)]
    expect_identical(fileNames(subs), fileNames(onDisk)[2])
    expect_true(all(fromFile(subs) == 1))
    expect_identical(pData(subs), droplevels(pData(onDisk)[2, , drop = FALSE]))
    expect_identical(experimentData(subs)@instrumentManufacturer,
                     experimentData(onDisk)@instrumentManufacturer[2])
    ## The same for MSnExp:
    subs <- inMem[c(1, 3, 5)]
    expect_identical(fileNames(subs), fileNames(inMem)[1])
    expect_true(all(fromFile(subs) == 1))
    expect_identical(pData(subs), droplevels(pData(inMem)[1, , drop = FALSE]))
    expect_identical(experimentData(subs)@instrumentManufacturer,
                     experimentData(inMem)@instrumentManufacturer[1])
    ## Extract spectra from the second file
    subs <- inMem[c(2, 4, 6)]
    expect_identical(fileNames(subs), fileNames(inMem)[2])
    expect_true(all(fromFile(subs) == 1))
    expect_identical(pData(subs), droplevels(pData(inMem)[2, , drop = FALSE]))
    expect_identical(experimentData(subs)@instrumentManufacturer,
                     experimentData(inMem)@instrumentManufacturer[2])
})
