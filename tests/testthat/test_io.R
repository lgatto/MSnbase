context("IO testing")

test_that("readMSData >1 files", {
    file1 <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
    tdir <- tempdir()
    file2 <- file.path(tdir, sub("iTRAQ", "iTRAQ2", basename(file1)))
    stopifnot(file.copy(file1, file2, overwrite = TRUE))
    raw <- readMSData(c(file1, file2))
    raw
})


test_that("Compatibility between writeMgfData and readMgfData", {
    data(itraqdata)
    tf <- tempfile()
    ## no spectra order issue here, as < 10 spectra
    ## see [read|write]MgfData for details
    d1 <- itraqdata[1:3]
    writeMgfData(d1, con = tf)
    d2 <- readMgfData(tf, verbose = FALSE)
    expect_true(all(sapply(1:3,
                           function(i) all.equal(as.data.frame(d1[[i]]),
                                                 as.data.frame(d1[[i]])))))
    expect_true(all.equal(precursorMz(d1), precursorMz(d1)))
    unlink(tf)
})

## Removing this test, as readMzXML has been defunct
## in MSnbase version 1.3.5
## test_that("Compare readMzXMLData and readMSData output", {
##   file <- dir(system.file(package="MSnbase", dir="extdata"),
##               full.names=TRUE,
##               pattern="mzXML$")
##   expect_warning(aa <- readMzXMLData(file, verbose=FALSE))
##   bb <- readMSData(file, verbose=FALSE)  
##   ## comparing data
##   expect_true(all.equal(header(aa), header(bb)))
##   expect_true(all.equal(intensity(aa), intensity(bb)))
##   expect_true(all.equal(mz(aa), mz(bb)))  
## })

test_that("Testing write.exprs and readMSnSet", {
    data(itraqdata)
    colchars <- c("ProteinAccession", "PeptideSequence",
                  "retention.time", "precursor.mz")
    tf <- tempfile()
    x <- quantify(itraqdata, reporters = iTRAQ4,
                  BPPARAM = SerialParam(),
                  method = "max", verbose = FALSE)
    write.exprs(x, file = tf)
    y <- readMSnSet(tf)
    expect_true(all.equal(exprs(x), exprs(y)))
    ## unlink(tf)
    write.exprs(x, fDataCols = colchars, file = tf)
    tmp <- read.table(tf)
    expect_true(all(dim(tmp) == c(nrow(x), ncol(x) + length(colchars))))
    expect_true(all(colnames(tmp) == c(sampleNames(x), colchars)))
    expect_true(all(rownames(tmp) == featureNames(x)))
    unlink(tf)
})


test_that("readMSnSet2: MSnSet from a data.frame", {
    k <- data.frame(A = 1:10, B = 10:1,
                    X1 = LETTERS[1:10], X2 = letters[1:10],
                    row.names = paste0("X", 1:10))
    x <- readMSnSet2(k, ecol = 1:2)
    expect_true(validObject(x))
    expect_identical(sampleNames(x), c("A", "B"))
    expect_identical(featureNames(x), rownames(k))
    expect_equivalent(exprs(x)[, 1], 1:10)
    expect_equivalent(exprs(x)[, 2], 10:1)
    expect_identical(fData(x), k[, 3:4])
    ## feature names as a column
    k$fn <- paste0("P", 1:10)
    x <- readMSnSet2(k, ecol = 1:2, fnames = "fn")
    expect_identical(featureNames(x), k$fn)
    rownames(k) <- k$fn
    expect_identical(fData(x), k[, 3:5])
    expect_error(readMSnSet2(k, ecol = 1:2, fnames = "fnames"))
    ## no feature names
    k <- data.frame(A = 1:10, B = 10:1,
                    X1 = LETTERS[1:10], X2 = letters[1:10])
    x <- readMSnSet2(k, ecol = 1:2)
    expect_identical(featureNames(x), as.character(1:10))
    x2 <- readMSnSet2(k, ecol = c("A", "B"))
    expect_identical(exprs(x), exprs(x2))
    expect_identical(fData(x), fData(x2))
    expect_equal(x, x2)
})


test_that("readMSnSet2: rownames and fnames", {
    f0 <- dir(system.file("extdata", package = "pRolocdata"),
              full.names = TRUE,
              pattern = "hyperLOPIT-SIData-ms3-rep12-intersect.csv")
    res1 <- readMSnSet2(f0, ecol = 8:27, fnames = 1, skip = 1)
    res2 <- readMSnSet2(f0, ecol = 8:27, rownames = 1, skip = 1)
    expect_warning(res3 <- readMSnSet2(f0, ecol = 8:27,
                                       rownames = 1, fnames = 1, skip = 1))
    expect_equal(res1, res2)
    expect_equal(res1, res3)
    f0 <- dir(system.file("extdata", package = "pRolocdata"),
              full.names = TRUE, pattern = "Dunkley2006")
    res1 <- readMSnSet2(f0, ecol = 5:20, fnames = 1)
    res2 <- readMSnSet2(f0, ecol = 5:20, fnames = "Protein.ID")
    res3 <- readMSnSet2(f0, ecol = 5:20, rownames = 1)
    res4 <- readMSnSet2(f0, ecol = 5:20, rownames = "Protein.ID")
    expect_equal(res1, res2)
    expect_equal(res3, res4)
    expect_equal(res1, res3)
})
