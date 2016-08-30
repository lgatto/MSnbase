context("Quantitation")

test_that("MS2 isobaric quantitation", {
    ## removeReporters
    expect_true(all.equal(removeReporters(itraqdata[[1]]),
                          itraqdata[[1]]))
    expect_true(all.equal(removeReporters(itraqdata,
                                          reporters = iTRAQ4,
                                          verbose = FALSE)[[1]],
                          removeReporters(itraqdata[[1]],
                                          reporters = iTRAQ4)))
    ## quantitation should be 0
    expect_true(all(quantify(removeReporters(itraqdata[[1]],
                                             reporters = iTRAQ4),
                             method = "max", reporters = iTRAQ4)[[1]] == 0))
    ## checking that quantification work for exp of length 1
    q1 <- quantify(itraqdata[1], reporters = iTRAQ4, method = "max",
                   BPPARAM = SerialParam())
    expect_true(inherits(q1, "MSnSet"))
})

    
test_that("Parallel quantification", {
    q1 <- quantify(itraqdata[1:5], reporters = iTRAQ4,
                   parallel = TRUE, verbose = FALSE)
    q2 <- quantify(itraqdata[1:5], reporters = iTRAQ4,
                   parallel = FALSE, verbose = FALSE)
    q1@processingData <- q2@processingData ## those are not expected to be equal
    expect_true(all.equal(q1, q2))
})


test_that("Counting and tic MSnSets", {
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    ## count
    .cnt <- MSnbase:::count_MSnSet(msx)
    n <- !is.na(fData(msx)$pepseq)
    m1 <- matrix(1, nrow = sum(n), ncol = 1)
    colnames(m1) <- "dummyiTRAQ.mzXML"
    rownames(m1) <- paste0("X", (1:length(msx))[n], ".1")
    expect_equal(exprs(.cnt), m1)
    ## tic
    .tic <- MSnbase:::tic_MSnSet(msx)
    mtic <- matrix(tic(msx), ncol = 1)
    colnames(mtic) <- "1"
    colnames(mtic) <- "dummyiTRAQ.mzXML"
    rownames(mtic) <- paste0("X", 1:length(msx), ".1")
    expect_equal(exprs(.tic), mtic)  
})

test_that("MS2 labelfree quantitation: SI", {
    ## prepare data
    f <- dir(system.file(package = "MSnbase",dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    fData(msx)$accession[3:4] <- "protein"
    fData(msx)$nprot[3:4] <- 1
    fData(msx)$pepseq[3:4] <- c("ABCDEFG", "1234567")
    fData(msx)$length[3:4] <- 100   
    msx <- MSnbase:::utils.removeNoIdAndMultipleAssignments(msx)
    ## SI
    si <- MSnbase:::SI(msx, "SI")
    m <- tic(msx)
    names(m) <- fData(msx)$accession
    expect_equal(exprs(si)["ECA0510", 1], as.numeric(m["ECA0510"]))
    expect_equal(exprs(si)["ECA1028", 1], as.numeric(m["ECA1028"]))
    expect_equal(exprs(si)["protein", 1], sum(m[2:3]))
    ## SIgi
    sigi <- MSnbase:::SI(msx, "SIgi")
    m <- c(m[1], protein = sum(m[2:3]), m[4])
    m <- m/sum(m)
    m <- m[order(names(m))]
    expect_equal(exprs(sigi)[, 1], m)
    ## SIn
    sin <- MSnbase:::SI(msx, "SIn")
    m <- m/fData(sin)$length
    expect_equal(exprs(sin)[, 1], m)
})

test_that("MS2 labelfree quantitation: SAF", {
    ## prepare data
    f <- dir(system.file(package = "MSnbase", dir = "extdata"),
             full.name = TRUE, pattern = "msx.rda")
    load(f) ## msx
    fData(msx)$accession[3:4] <- "protein"
    fData(msx)$nprot[3:4] <- 1
    fData(msx)$pepseq[3:4] <- c("ABCDEFG", "1234567")
    fData(msx)$length[3:4] <- 100   
    msx <- MSnbase:::utils.removeNoIdAndMultipleAssignments(msx)
    ## SAF
    saf <- MSnbase:::SAF(msx, method = "SAF")
    m <- rep(1, length(msx))
    m <- tapply(m, fData(msx)$accession, sum, simplify = FALSE)
    m <- unlist(m)
    expect_equal(featureNames(saf), names(m))
    m <- m/fData(saf)$length
    expect_equal(exprs(saf)[, 1], m)
    ## NSAF
    nsaf <- MSnbase:::SAF(msx, method = "NSAF")
    m <- m/sum(m)
    expect_equal(exprs(nsaf)[, 1], m)
})


test_that("quantify_OnDiskMSnExp_max (fastquant_max)", {
    f <- dir(system.file(package = "MSnbase",dir = "extdata"),
             full.name = TRUE, pattern = "dummyiTRAQ.mzXML")
    x1 <- readMSData(f, msLevel = 2, verbose = FALSE)
    x2 <- readMSData2(f, msLevel = 2, verbose = FALSE)
    e1 <- quantify(x1, method = "max", reporters = iTRAQ4, verbose = FALSE)
    e2 <- quantify(x2, method = "max", reporters = iTRAQ4, verbose = FALSE)
    expect_identical(exprs(e1), exprs(e2))
})
