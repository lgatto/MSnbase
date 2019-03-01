context("BackendMemory class")

test_that("constructor", {
    b <- BackendMemory()
    expect_s4_class(b, "BackendMemory")
    expect_length(b@spectra, 0)
})

test_that("validity", {
    b <- BackendMemory()
    b@spectra <- c("F1.S1"=new("Spectrum2"), "F2.S2"=new("Spectrum2"))
    b@files <- "foo"
    expect_error(validObject(b), "counters")
    b@modCount <- 1L
    expect_true(validObject(b))

    names(b@spectra)[2] <- "F1.S1"
    expect_error(validObject(b), "Duplicated spectra names")

    names(b@spectra)[2] <- "F1.S2"
    b@spectra[[1]]@peaksCount <- 1L
    expect_error(validObject(b), "Peaks count does not match")
})

test_that(".valid.BackendMemory.spectra.names", {
    expect_null(.valid.BackendMemory.spectra.names(list()))
    expect_null(.valid.BackendMemory.spectra.names(list(F1.S1=1, F1.S2=2)))
    expect_match(.valid.BackendMemory.spectra.names(list(1, 2)), " NULL")
    expect_match(.valid.BackendMemory.spectra.names(list(a=1, 2)), " missing")
    l <- list(1, 2)
    names(l)[1] <- "a"
    expect_match(.valid.BackendMemory.spectra.names(l), " NA")
    expect_match(.valid.BackendMemory.spectra.names(list(a=1, a=2)),
                 "Duplicated")
    ## expect_match(.valid.BackendMemory.spectra.names(list(a=1, b=2)),
    ##              "format")
})

test_that("backendSubset,BackendMemory works", {
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    tmp <- readMSnExperiment(f, backend = BackendMemory())
    be <- tmp@backend
    spd <- tmp@spectraData
    sps <- spectrapply(tmp)
    ## Subset to data from the second file.
    be_2 <- backendSubset(be, spd[spd$fileIdx == 2, ])
    expect_equal(be_2@files, be@files[2])
    ## fromFile has to be 1 for all spectra
    expect_true(all(vapply(be_2@spectra, fromFile, integer(1)) == 1))
    expect_equal(lapply(be_2@spectra, intensity),
                 lapply(sps[spd$fileIdx == 2], intensity))
    ## Subset to some specific spectra.
    idx <- c(200, 201, 3, 5, 6)
    be_3 <- backendSubset(be, spd[idx, ])
    expect_equal(be_3@files, be@files[2:1])
    expect_equal(unname(vapply(be_3@spectra, fromFile, integer(1))),
                 c(1, 1, 2, 2, 2))
    expect_equal(lapply(be_3@spectra, intensity), lapply(sps[idx], intensity))
})

test_that("backendSplitByFile,BackendMemory works", {
    b <- BackendMemory()
    b@files <- c("a", "b", "c")
    s <- c(F1.S1=new("Spectrum2", mz=1:2, intensity=1:2, fromFile = 1L),
           F2.S2=new("Spectrum2", mz=3:4, intensity=3:4, fromFile = 2L),
           F3.S3=new("Spectrum2", mz=5:6, intensity=5:6, fromFile = 3L),
           F1.S2=new("Spectrum2", mz=7:8, intensity=7:8, fromFile = 1L))
    b@spectra <- s[1:3]
    b@modCount <- rep(0L, 3L)
    spd <- DataFrame(fileIdx = c(3, 1, 2))
    rownames(spd) <- names(s)[c(3, 1, 2)]
    res <- backendSplitByFile(b, spd)
    bl <- BackendMemory()
    bl@files <- "a"
    bl@spectra <- s[1]
    bl@modCount <- 0L
    l <- list("1"=bl, "2"=bl, "3"=bl)
    l[[2]]@files <- "b"
    l[[2]]@spectra <- s[2]
    l[[2]]@spectra[[1]]@fromFile <- 1L
    l[[3]]@files <- "c"
    l[[3]]@spectra <- s[3]
    l[[3]]@spectra[[1]]@fromFile <- 1L
    expect_equal(backendSplitByFile(b, spd), l)
    r <- b
    r@files[1] <- "d"
    r@spectra[1] <- s[4]
    r@modCount[1L] <- 1L
    l[[1]]@files <- "d"
    l[[1]]@spectra <- s[4]
    l[[1]]@modCount <- 1L
    backendSplitByFile(b, spd) <- l
    expect_equal(b, r)
})

test_that("backendInitialize", {
    spd <- DataFrame(fileIdx=c(1, 1, 2), spIdx=1:3,
                     row.names=c("F1.S1", "F1.S2", "F2.S3"))
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    b <- backendInitialize(BackendMemory(), files=f, spectraData=spd)
    expect_length(b@spectra, 3)
    expect_equal(names(b@spectra), rownames(spd))
})

test_that("backendReadSpectra/backendWriteSpectra", {
    spd <- DataFrame(fileIdx=c(1, 1, 2), spIdx=1:3,
                     row.names=c("F1.S1", "F1.S2", "F2.S3"))
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    b <- backendInitialize(BackendMemory(), files=f, spectraData=spd)
    s <- c(F1.S1=new("Spectrum2", mz=1:2, intensity=1:2, fromFile = 1L),
           F1.S2=new("Spectrum2", mz=3:4, intensity=3:4, fromFile = 1L),
           F2.S3=new("Spectrum2", mz=5:6, intensity=5:6, fromFile = 2L))
    b@spectra[] <- s
    expect_equal(backendReadSpectra(b, spd[1:2,]), s[1:2])
    expect_equal(backendReadSpectra(b, spd[3,]), s[3])

    r <- b
    r@spectra[] <- s[c(1, 2, 2)]
    r@modCount <- c(1L, 0L)
    expect_equal(backendWriteSpectra(b, s[2], spd[3,]), r)
    r@spectra[] <- s[c(2, 1, 3)]
    r@modCount <- c(1L, 0L)
    expect_equal(backendWriteSpectra(b, s[2:1], spd[1:2,]), r)

    r@spectra[] <- s
    r@modCount <- c(1L, 1L)
    expect_equal(backendWriteSpectra(b, s, spd), r)
})

test_that("backendUpdateMetadata,BackendMemory works", {
    spd <- DataFrame(fileIdx=c(1, 1, 2), spIdx=1:3,
                     row.names=c("F1.S1", "F1.S2", "F2.S3"))
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    b <- backendInitialize(BackendMemory(), files=f, spectraData=spd)
    s <- c(F1.S1=new("Spectrum2", mz=1:2, intensity=1:2),
           F1.S2=new("Spectrum2", mz=3:4, intensity=3:4),
           F2.S3=new("Spectrum2", mz=5:6, intensity=5:6))
    b@spectra[] <- s
    spd$msLevel <- c(2, 3, 4)
    spd$polarity <- c(-1, -1, 1)
    res <- backendUpdateMetadata(b, spd)
    expect_equal(polarity(res@spectra[[1]]), -1)
    expect_equal(polarity(res@spectra[[2]]), -1)
    expect_equal(polarity(res@spectra[[3]]), 1)
    expect_equal(msLevel(res@spectra[[1]]), 2)
    expect_equal(msLevel(res@spectra[[2]]), 3)
    expect_equal(msLevel(res@spectra[[3]]), 4)
})
