context("utils")

test_that("getColsFromPattern", {
  m <- matrix(1:12, nrow=4)
  expect_error(MSnbase:::getColsFromPattern(1:10, "101"),
               "must be a matrix")
  expect_error(MSnbase:::getColsFromPattern(m),
               ".*pattern.* must not be missing")
  expect_error(MSnbase:::getColsFromPattern(m, "1011"),
               "must be equal to the number of columns")
  expect_equal(MSnbase:::getColsFromPattern(m, "111"), rep(TRUE, 3))
  expect_equal(MSnbase:::getColsFromPattern(m, "000"), rep(FALSE, 3))
  expect_equal(MSnbase:::getColsFromPattern(m, "101"), c(TRUE, FALSE, TRUE))
  expect_equal(MSnbase:::getColsFromPattern(m, "010"), c(FALSE, TRUE, FALSE))
})

test_that("getRowsFromPattern", {
  m <- matrix(1:12, nrow=4)
  m[c(2, 12)] <- NA_real_
  expect_error(MSnbase:::getRowsFromPattern(1:10, "101"),
               "must be a matrix")
  expect_error(MSnbase:::getRowsFromPattern(m),
               ".*pattern.* must not be missing")
  expect_error(MSnbase:::getRowsFromPattern(m, "1011"),
               "must be equal to the number of columns")
  expect_equal(MSnbase:::getRowsFromPattern(m, "000"), rep(TRUE, 4))
  expect_equal(MSnbase:::getRowsFromPattern(m, "111"),
               c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::getRowsFromPattern(m, "101"),
               c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::getRowsFromPattern(m, "010"),
               rep(TRUE, 4))
  expect_equal(MSnbase:::getRowsFromPattern(m, "110"),
               c(TRUE, FALSE, TRUE, TRUE))
})

test_that(".filterNA", {
  m <- matrix(1:12, nrow=4)
  m[c(2, 12)] <- NA_real_
  expect_error(MSnbase:::.filterNA(1:10, 1, "110"), "must be a matrix")
  expect_error(MSnbase:::.filterNA(m, pNA="110"), "must be numeric")
  expect_error(MSnbase:::.filterNA(m, 1:10), "must be of length one")
  expect_equal(MSnbase:::.filterNA(m), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::.filterNA(m, pNA=1), c(TRUE, TRUE, TRUE, TRUE))
  expect_equal(MSnbase:::.filterNA(m, pNA=0.4), c(TRUE, TRUE, TRUE, TRUE))
  expect_equal(MSnbase:::.filterNA(m, pNA=0.3), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::.filterNA(m, pattern="101"), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("vec2ssv & ssv2vec", {
  numbers <- 1:3
  string1 <- "1;2;3"
  string2 <- "1,2,3"

  expect_equal(MSnbase:::utils.vec2ssv(numbers), string1)
  expect_equal(MSnbase:::utils.vec2ssv(numbers, sep=","), string2)

  expect_equal(as.numeric(MSnbase:::utils.ssv2vec(string1)), numbers)
  expect_equal(as.numeric(MSnbase:::utils.ssv2vec(string2, sep=",")), numbers)
})

test_that("list2ssv & ssv2list", {
    l <- list(a=1:3, b=4:6)
    string1 <- c("1;2;3;4;5;6")
    string2 <- c("1,2,3,4,5,6")
    strings1 <- c(a="1;2;3", b="4;5;6")
    strings2 <- c(a="1,2,3", b="4,5,6")

    expect_equal(MSnbase:::utils.list2ssv(l), strings1)
    expect_equal(MSnbase:::utils.list2ssv(l, sep=","), strings2)

    expect_equal(lapply(MSnbase:::utils.ssv2list(strings1), as.numeric), l)
    expect_equal(lapply(MSnbase:::utils.ssv2list(strings2, sep=","), as.numeric), l)

    expect_equal(MSnbase:::utils.vec2ssv(unlist(l)), string1)
    expect_equal(MSnbase:::utils.vec2ssv(unlist(l), sep=","), string2)
})

test_that("leftJoin", {
    x <- data.frame(id=1:7,
                    fn=LETTERS[1:7],
                    add=1:7,
                    useless1=1,
                    stringsAsFactors=FALSE)
    rownames(x) <- paste0("R", seq(nrow(x)))
    y1 <- data.frame(id=4:2,
                     fn=LETTERS[4:2],
                     foobar=letters[4:2],
                     useless2=2,
                     stringsAsFactors=FALSE)
    y2 <- data.frame(id=6:7,
                     fn=LETTERS[6:7],
                     foobar=letters[6:7],
                     useless3=3,
                     stringsAsFactors=FALSE)
    z1 <- data.frame(id=1:7,
                     fn=LETTERS[1:7],
                     add=1:7,
                     foobar=c(NA, letters[2:4], rep(NA, 3)),
                     stringsAsFactors=FALSE)
    z2 <- data.frame(id=1:7,
                     fn=LETTERS[1:7],
                     add=1:7,
                     foobar=c(NA, letters[2:4], NA, letters[6:7]),
                     stringsAsFactors=FALSE)
    rownames(z1) <- rownames(z2) <- paste0("R", seq(nrow(x)))
    ## first run
    expect_equal(MSnbase:::utils.leftJoin(x, y1, by=c("id", "fn"),
                                          exclude=c("useless1", "useless2")), z1)
    ## second run (on the results of the first run)
    expect_equal(MSnbase:::utils.leftJoin(z1, y2, by=c("id", "fn"),
                                          exclude=c("useless3")), z2)
})

test_that("mergeSpectraAndIdentificationData", {
    ## pseudo fData(MSnSet) output
    fd <- data.frame(spectrum = 1:4,
                     file = c(1, 2, 1, 1),
                     acquisition.number = 5:8,
                     row.names = paste0("R", 1:4),
                     stringsAsFactors = FALSE)
    ## pseudo mzID output
    id1 <- data.frame(acquisitionnum = c(5, 5, 5, 8),
                      file = 1,
                      spectrumFile = "foobar1.mzML",
                      rank = c(2, 3, 1, 1),
                      accession = paste0("P", 1:4),
                      description = paste0("D", 1:4),
                      spectrumid = paste0("id", 1:4),
                      pepseq = LETTERS[1:4],
                      stringsAsFactors = FALSE)
    id2 <- data.frame(acquisitionnum = 6,
                      file = 2,
                      spectrumFile = "foobar2.mzML",
                      rank = 1,
                      accession = "P9",
                      description = "D9",
                      spectrumid = "id9",
                      pepseq = "F",
                      stringsAsFactors = FALSE)
    ## results
    rfd1 <- data.frame(spectrum = 1:4,
                       file = c(1, 2, 1, 1),
                       acquisition.number = 5:8,
                       rank = c(1, NA, NA, 1),
                       accession = c("P3;P1;P2", NA, NA, "P4"),
                       description = c("D3;D1;D2", NA, NA, "D4"),
                       pepseq = c("C", NA, NA, "D"),
                       nprot = c(3, NA, NA, 1),
                       npep.prot = c(1, NA, NA, 1),
                       npsm.prot = c(1, NA, NA, 1),
                       npsm.pep = c(1, NA, NA, 1),
                       row.names = paste0("R", 1:4),
                       stringsAsFactors = FALSE)
    rfd2 <- data.frame(spectrum = 1:4,
                       file = c(1, 2, 1, 1),
                       acquisition.number = 5:8,
                       rank = c(1, 1, NA, 1),
                       accession = c("P3;P1;P2", "P9", NA, "P4"),
                       description = c("D3;D1;D2", "D9", NA, "D4"),
                       pepseq = c("C", "F", NA, "D"),
                       nprot = c(3, 1, NA, 1),
                       npep.prot = c(1, 1, NA, 1),
                       npsm.prot = c(1, 1, NA, 1),
                       npsm.pep = c(1, 1, NA, 1),
                       row.names = paste0("R", 1:4),
                       stringsAsFactors = FALSE)
    ## errors
    expect_error(MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id1,
                                                                   fcol = c("file", "acquisition_number"),
                                                                   icol = c("file", "acquisitionnum")))
    expect_error(MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id1,
                                                                   fcol = c("file", "acquisition.number"),
                                                                   icol = c("file", "acquisition_num")))
    expect_error(MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id1,
                                                                   fcol = c("file", "acquisition.number"),
                                                                   icol = c("file", "acquisitionnum", "rank")))
    ## first run
    ans1 <-
        MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id1,
                                                          fcol = c("file", "acquisition.number"),
                                                          icol = c("file", "acquisitionnum"),
                                                          acc = "accession",
                                                          desc = "description", pepseq = "pepseq")
    expect_equal(ans1, rfd1)    
    ## second run
    ans2 <-
        MSnbase:::utils.mergeSpectraAndIdentificationData(rfd1, id2,
                                                          fcol = c("file", "acquisition.number"),
                                                          icol = c("file", "acquisitionnum"),
                                                          acc = "accession",
                                                          desc = "description", pepseq = "pepseq")
    expect_equal(ans2, rfd2)    
})

test_that("utils.idSummary", {
    ## pseudo fData(MSnSet) output
    fd <- data.frame(spectrumFile = c(1, 2, 1, 1, 1, 3),
                     idFile = c(4, 5, 4, 4, NA, NA))
    ## results
    rdf <- data.frame(spectrumFile = c(1, 2, 3),
                      idFile = c(4, 5, NA),
                      coverage = c(0.75, 1, 0))
    expect_error(MSnbase:::utils.idSummary(data.frame(file = 1:3, foobar = 1:3)),
                 "No quantification/identification data found")
    expect_equal(MSnbase:::utils.idSummary(fd), rdf)
})

test_that("formatRt", {
    tc <- c("1:1", "25:24")
    tn <- c(61, 25 * 60 + 24)
    expect_equal(tc, formatRt(tn))
    expect_equal(tn, formatRt(tc))
    expect_true(is.na(formatRt("")))
    expect_warning(is.na(formatRt("aaa")))
    expect_warning(is.na(formatRt(TRUE)))
})

test_that("formatFileSpectrumNames", {
    expect_error(MSnbase:::formatFileSpectrumNames(1:10, 1:3),
                 "has to be one or equal")
    expect_equal(MSnbase:::formatFileSpectrumNames(1, 1), "F1.S1")
    expect_equal(MSnbase:::formatFileSpectrumNames(1, 1:10),
                 sprintf("F1.S%02d", 1:10))
    expect_equal(MSnbase:::formatFileSpectrumNames(1:10, 1:10),
                 sprintf("F%02d.S%02d", 1:10, 1:10))
    expect_equal(MSnbase:::formatFileSpectrumNames(1:10, 1:10, nSpectra=1000),
                 sprintf("F%02d.S%04d", 1:10, 1:10))
    expect_equal(MSnbase:::formatFileSpectrumNames(1:10, 1:10,
                                                   nSpectra=1000, nFiles=1000),
                 sprintf("F%04d.S%04d", 1:10, 1:10))
})

test_that("rowmean", {
    m <- matrix(1:10, ncol=2)
    mna <- m
    mna[c(2, 8)] <- NA
    group <- c("B", "B", "A", "A", "B")
    r <- matrix(c(8/3, 3.5, 23/3, 8.5), ncol=2, dimnames=list(c("B", "A"), c()))
    rna <- r
    rna[c(1, 4)] <- NA
    rnarm <- matrix(c(3, 3.5, 23/3, 9), ncol=2, dimnames=list(c("B", "A"), c()))
    expect_error(MSnbase:::rowmean(m, group=1:2))
    expect_equal(MSnbase:::rowmean(m, group=group), r)
    expect_equal(MSnbase:::rowmean(m, group=group, reorder=TRUE), r[2:1,])
    expect_equal(MSnbase:::rowmean(mna, group=group), rna)
    expect_equal(MSnbase:::rowmean(mna, group=group, reorder=TRUE), rna[2:1,])
    expect_equal(MSnbase:::rowmean(mna, group=group, na.rm=TRUE), rnarm)
})

test_that("rowsd", {
    m <- matrix(1:10, ncol=2)
    mna <- m
    mna[c(2, 8)] <- NA
    group <- c("B", "B", "A", "A", "B")
    r <- matrix(c(sd(m[c(1, 2, 5), 1]), sd(m[3:4, 1]),
                  sd(m[c(1, 2, 5), 2]), sd(m[3:4, 2])),
                  ncol=2, dimnames=list(c("B", "A"), c()))
    rna <- r
    rna[c(1, 4)] <- NA
    rnarm <- matrix(c(sd(mna[c(1, 2, 5), 1], na.rm=TRUE), r[2],
                      r[3], sd(mna[3:4, 2], na.rm=TRUE)),
          ncol=2, dimnames=list(c("B", "A"), c()))
    expect_error(MSnbase:::rowsd(m, group=1:2))
    expect_equal(MSnbase:::rowsd(m, group=group), r)
    expect_equal(MSnbase:::rowsd(m, group=group, reorder=TRUE), r[2:1,])
    expect_equal(MSnbase:::rowsd(mna, group=group), rna)
    expect_equal(MSnbase:::rowsd(mna, group=group, reorder=TRUE), rna[2:1,])
    expect_equal(MSnbase:::rowsd(mna, group=group, na.rm=TRUE), rnarm)
})

test_that("get.amino.acids", {
    aa <- get.amino.acids()
    cn <- c("AA", "ResidueMass", "Abbrev3", "ImmoniumIonMass", "Name",
            "Hydrophobicity", "Hydrophilicity", "SideChainMass",
            "pK1", "pK2", "pI") ## from man
    expect_identical(colnames(aa), cn)
})

test_that("remove precursor MZ", {
    data(itraqdata)
    sp <- itraqdata[[1]]
    r1 <- MSnbase:::utils.removePrecMz(sp)
    spl <- list(mz = mz(sp),
                int = intensity(sp))
    r2 <- MSnbase:::utils.removePrecMz_list(spl,
                                            precursorMz(sp))
    expect_identical(mz(r1), r2$mz)

    r1 <- MSnbase:::utils.removePrecMz(sp, width = 1)
    spl <- list(mz = mz(sp),
                int = intensity(sp))
    r2 <- MSnbase:::utils.removePrecMz_list(spl,
                                            precursorMz(sp),
                                            width = 1)
    expect_identical(mz(r1), r2$mz)

    spl <- list(mz = 1:10, int = 1:10)
    spl <- MSnbase:::utils.removePrecMz_list(spl, 5, 2)
    expect_identical(spl$mz, 1:10)
    expect_identical(spl$int, c(1:3, rep(0, 3), 7:10))
})


test_that("clean utils", {
    x <- list(c(0, 0, 0, 1, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
              c(1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0),
              c(1, NA, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, NA, 1, 0, 0, 0))
    a <- list(c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
                FALSE),
              c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
                FALSE, FALSE), c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
                TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
              c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
                TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
    b <- list(c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
    d <- list(c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                FALSE, FALSE),
              c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
              c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
    for (i in seq(along=x)) {
      expect_identical(MSnbase:::utils.clean(x[[i]], all=TRUE), a[[i]],
                       label=paste0("all=TRUE, i=", i))
      expect_identical(MSnbase:::utils.clean(x[[i]], all=FALSE), b[[i]],
                       label=paste0("all=FALSE, i=", i))
      expect_identical(MSnbase:::utils.clean(x[[i]], na.rm=TRUE), d[[i]],
                       label=paste0("na.rm=TRUE, i=", i))
    }
})

test_that("enableNeighbours utils", {
    expect_error(MSnbase:::utils.enableNeighbours(1:10))
    expect_error(MSnbase:::utils.enableNeighbours(LETTERS[1:10]))
    expect_equal(MSnbase:::utils.enableNeighbours(c(FALSE, TRUE, FALSE)),
                 rep(TRUE, 3))
    expect_equal(MSnbase:::utils.enableNeighbours(c(FALSE, TRUE, FALSE, FALSE)),
                 c(rep(TRUE, 3), FALSE))
})

test_that("getBpParam", {
    ## Testing the global MSnbase option PARALLEL_THRESH
    orig_val <- options()$MSnbase$PARALLEL_THRESH
    suppressWarnings(
        onDisk <- readMSData2(files = system.file("microtofq/MM14.mzML",
                                                  package = "msdata"),
                              verbose = FALSE)
    )
    gotParam <- MSnbase:::getBpParam(onDisk)
    expect_true(is(gotParam, "SerialParam"))
    options(MSnbase=list(PARALLEL_THRESH = 10))
    gotParam <- MSnbase:::getBpParam(onDisk)
    expect_true(is(gotParam, class(bpparam())))
    options(MSnbase=list(PARALLEL_THRESH = orig_val))
})

test_that("Get first MS level", {
    MSnbase::setMSnbaseVerbose(FALSE)
    f <- msdata::proteomics(full.names = TRUE, pattern = "MS3TMT10")
    x <- readMSData(f, msLevel. = 2L)
    y <- readMSData2(f, msLevel. = 2L)
    ## in memory
    tx1 <- system.time(x1 <- msLevel(x)[1])[["elapsed"]]
    tx2 <- system.time(x2 <- MSnbase:::.firstMsLevel(x))[["elapsed"]]
    expect_equivalent(x1, x2)
    ## expect_lt(tx2, tx1)
    ## on disk - no timing here, as msLevel is extracted from fData
    ## and hence already fast; same implementation in msLevel and
    ## .firstMsLevel
    y1 <- msLevel(y)[1]
    y2 <- MSnbase:::.firstMsLevel(y)
    expect_equivalent(y1, y2)
})

test_that(".rowMaxs", {
  m <- matrix(1:30, nrow=10)
  m[c(2, 5:6, 10, 12:14, 21:25)] <- NA
  expect_error(MSnbase:::.rowMaxs(1:10))
  expect_equal(MSnbase:::.rowMaxs(m), apply(m, 1, max))
  expect_equal(MSnbase:::.rowMaxs(m, na.rm=TRUE),
               suppressWarnings(apply(m, 1, max, na.rm=TRUE)))
})

test_that(".summariseRows", {
  m <- matrix(1:30, nrow=10)
  m[seq(2, 30, by=3)] <- NA
  expect_error(MSnbase:::.summariseRows(1:10, fun=sum))
  expect_error(MSnbase:::.summariseRows(m, fun=1:10))
  expect_equal(MSnbase:::.summariseRows(m, fun=sum), rep.int(NA_real_, 10))
  expect_equal(MSnbase:::.summariseRows(m, fun=sum, na.rm=TRUE), rowSums(m, na.rm=TRUE))
  expect_equal(MSnbase:::.summariseRows(m, fun=mean, na.rm=TRUE), rowMeans(m, na.rm=TRUE))
  expect_equal(MSnbase:::.summariseRows(m, fun=max, na.rm=TRUE), c(21, 22, 13, 24, 25, 16, 27, 28, 19, 30))
})

test_that(".topIdx", {
  m <- matrix(1:30, nrow=10)
  g <- rep_len(LETTERS[1:3], 10)
  expect_error(MSnbase:::.topIdx(m, groupBy=g, fun=sum)) # n missing
  expect_error(MSnbase:::.topIdx(1:10, groupBy=g, fun=sum, n=3))
  expect_error(MSnbase:::.topIdx(m, groupBy=g, fun=1:10, n=3))
  expect_error(MSnbase:::.topIdx(m, groupBy=g, fun=sum, n=-1),
               ".*n.* has to be greater or equal than 1.")
  expect_error(MSnbase:::.topIdx(m, groupBy=1:3, fun=sum, n=3),
               ".*nrow.*x.* and .*length.*groupBy.* have to be equal.")
  expect_equal(MSnbase:::.topIdx(m, groupBy=g, fun=sum, n=3),
               c(10, 7, 4, 8, 5, 2, 9, 6, 3))
  expect_equal(MSnbase:::.topIdx(m, groupBy=g, fun=sum, n=2),
               c(10, 7, 8, 5, 9, 6))
})

test_that(".mzRBackend works", {
    res <- MSnbase:::.mzRBackend("test.mzml")
    expect_equal(res, "pwiz")
    res <- MSnbase:::.mzRBackend("test.mzml.gz")
    expect_equal(res, "pwiz")
    res <- MSnbase:::.mzRBackend("test.mzML")
    expect_equal(res, "pwiz")
    res <- MSnbase:::.mzRBackend("test.mzML.bz2")
    expect_equal(res, "pwiz")
    res <- MSnbase:::.mzRBackend("test.mzXML")
    expect_equal(res, "pwiz")

    res <- MSnbase:::.mzRBackend("test.mzdata")
    expect_equal(res, "Ramp")
    res <- MSnbase:::.mzRBackend("test.mzdata.gz")
    expect_equal(res, "Ramp")
    
    res <- MSnbase:::.mzRBackend("test.cdf")
    expect_equal(res, "netCDF")
    res <- MSnbase:::.mzRBackend("test.cdf.gz")
    expect_equal(res, "netCDF")

    expect_error(MSnbase:::.mzRBackend("unsupported.txt"))
    expect_error(MSnbase:::.mzRBackend())
    expect_error(MSnbase:::.mzRBackend(""))
    expect_error(MSnbase:::.mzRBackend(c("a.mzML", "b.mzML")))
})

test_that(".openMSfile works", {
    file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    res <- MSnbase:::.openMSfile(file)
    expect_true(is(res, "mzRpwiz"))
    close(res)
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                        package = "msdata")
    res <- MSnbase:::.openMSfile(file)
    expect_true(is(res, "mzRpwiz"))
    close(res)
    file <- system.file("microtofq", "MM14.mzdata", package = "msdata")
    res <- MSnbase:::.openMSfile(file)
    expect_true(is(res, "mzRramp"))
    close(res)

    ## Errors
    expect_error(MSnbase:::.openMSfile(c("a", "b")))
    file <- c(system.file("microtofq", "MM14.mzML", package = "msdata"),
              system.file("microtofq", "MM14.mzdata", package = "msdata"))
    expect_error(MSnbase:::.openMSfile(file))
})

test_that("camel case", {
    k0 <- c("aa.bb", "cc.dd")
    expect_identical(makeCamelCase(k0), c("aaBb", "ccDd"))
    expect_identical(makeCamelCase(k0, prefix = "x"),
                     c("xAaBb", "xCcDd"))
})

test_that("factorsAsStrings",{
    data(iris)
    expect_true(is.factor(iris$Species))
    iris2 <- factorsAsStrings(iris)
    expect_true(is.character(iris2$Species))
    expect_identical(dim(iris), dim(iris2))
    expect_identical(iris[, -5], iris2[, -5])
})
