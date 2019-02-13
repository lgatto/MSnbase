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
    be_2 <- MSnbase:::backendSubset(be, spd[spd$fileIdx == 2, ])
    expect_equal(be_2@files, be@files[2])
    ## fromFile has to be 1 for all spectra
    expect_true(all(vapply(be_2@spectra, fromFile, integer(1)) == 1))
    expect_equal(lapply(be_2@spectra, intensity),
                 lapply(sps[spd$fileIdx == 2], intensity))
    ## Subset to some specific spectra.
    idx <- c(200, 201, 3, 5, 6)
    be_3 <- MSnbase:::backendSubset(be, spd[idx, ])
    expect_equal(be_3@files, be@files[2:1])
    expect_equal(unname(vapply(be_3@spectra, fromFile, integer(1))),
                 c(1, 1, 2, 2, 2))
    expect_equal(lapply(be_3@spectra, intensity), lapply(sps[idx], intensity))
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

test_that("backendImportData", {
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    nn <- c(4, 3)
    n <- sum(nn)
    spd <-  DataFrame(
        fileIdx=rep(1:2, nn),
        spIdx=c(seq_len(nn[1]), seq_len(nn[2])),
        smoothed=FALSE,
        row.names=paste0(rep(c("F1", "F2"), nn),
                         ".S", c(seq_len(nn[1]), seq_len(nn[2])))
    )
    hdr <- pks <- vector(mode="list", length=length(f))

    for (i in seq(along=f)) {
        fh <- mzR::openMSfile(f[i])
        j <- seq_len(nn[i])
        hdr[[i]] <- header(fh, j)
        pks[[i]] <- peaks(fh, j)
        close(fh)
    }
    hdr <- do.call(rbind, hdr)
    pks <- unlist(pks, recursive=FALSE)
    spd <- cbind(spd, hdr)

    b <- backendInitialize(BackendMemory(), files=f, spectraData=spd)
    b <- backendImportData(b, spectraData=spd, BPPARAM=SerialParam())
    expect_length(b@spectra, n)
    expect_equal(names(b@spectra), rownames(spd))
    expect_equal(lapply(b@spectra, mz),
                 setNames(lapply(pks, function(p)p[,1]), rownames(spd)))
    expect_equal(lapply(b@spectra, intensity),
                 setNames(lapply(pks, function(p)p[,2]), rownames(spd)))
})

test_that("backendReadSpectra/backendWriteSpectra", {
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
    expect_equal(backendReadSpectra(b, spd[1:2,]), s[1:2])
    expect_equal(backendReadSpectra(b, spd[3,]), s[3])

    r <- b
    r@spectra[] <- s[c(1, 2, 2)]
    expect_equal(backendWriteSpectra(b, s[2], spd[3,]), r)
    r@spectra[] <- s[c(2, 1, 3)]
    expect_equal(backendWriteSpectra(b, s[2:1], spd[1:2,]), r)
})

test_that("backendUpdateMetadata,BackendMemory works", {
    spd <- DataFrame(fileIdx=c(1, 1, 2), spIdx=1:3,
                     row.names=c("F1.S1", "F1.S2", "F2.S3"))
    f <- system.file(
        file.path("microtofq", c("MM8.mzML", "MM14.mzML")),
        package="msdata"
    )
    b <- MSnbase:::backendInitialize(BackendMemory(), files=f, spectraData=spd)
    s <- c(F1.S1=new("Spectrum2", mz=1:2, intensity=1:2),
           F1.S2=new("Spectrum2", mz=3:4, intensity=3:4),
           F2.S3=new("Spectrum2", mz=5:6, intensity=5:6))
    b@spectra[] <- s
    spd$msLevel <- c(2, 3, 4)
    spd$polarity <- c(-1, -1, 1)
    res <- MSnbase:::backendUpdateMetadata(b, spd)
    expect_equal(polarity(res@spectra[[1]]), -1)
    expect_equal(polarity(res@spectra[[2]]), -1)
    expect_equal(polarity(res@spectra[[3]]), 1)
    expect_equal(msLevel(res@spectra[[1]]), 2)
    expect_equal(msLevel(res@spectra[[2]]), 3)
    expect_equal(msLevel(res@spectra[[3]]), 4)
})
