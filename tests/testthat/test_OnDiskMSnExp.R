context("OnDiskMSnExp class")

.getMzMLFiles <- function(){
    ## Return the mzML files, the ones from the XXX package, or if run
    ## locally, some of my test files.
    HOST <- unlist(strsplit(system("hostname", intern=TRUE), split=".",
                            perl=FALSE, fixed=TRUE))[1]
    if(HOST == "macbookjo"){
        mzfiles <- dir("/Users/jo/R-workspaces/EURAC/2016/2016-04-21-PolarMetabolom/data/mzML/",
                       pattern="POS_C_O", full.names=TRUE)
    }else{
        require(msdata)
        mzfiles <- c(system.file("microtofq/MM14.mzML", package="msdata"),
                     system.file("microtofq/MM8.mzML", package="msdata"))
    }
    return(mzfiles)
}
mzf <- .getMzMLFiles()[1:2]
mse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="ram")
odmse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="disk")



############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
    expect_identical(as.character(class(mse)), "MSnExp")
    expect_identical(as.character(class(odmse)), "OnDiskMSnExp")
})

test_that("read and validate OnDiskMSnExp data", {
    ## Testing everything in one test fun, so we don't have to re-read the data
    ## again.
    ## Reading the stuff in memory.
    ## Read as an OnDiskMSnExp

    ## Check if we get the same data! MSnExp will retrieve that from the spectra,
    ## OnDiskMSnExp from featureData.
    ## fromFile
    expect_identical(fromFile(mse), fromFile(odmse))

    ## msLevel
    expect_identical(msLevel(mse), msLevel(odmse))

    ## header
    hd <- header(mse)
    odhd <- header(odmse)
    commonCols <- intersect(colnames(hd), colnames(odhd))
    expect_identical(hd[, commonCols], odhd[, commonCols])

    ## length
    expect_identical(length(mse), length(odmse))

    ## acquisitionNum
    expect_identical(acquisitionNum(mse), acquisitionNum(odmse))

    ## scanIndex; point is that scanIndex on an MSnExp will return 0, as it is not
    ## set in the spectra (same as acquisitionNum?)
    expect_identical(acquisitionNum(mse), scanIndex(odmse))

    ## centroided.
    expect_identical(centroided(mse), centroided(odmse))
    ## Setting stuff
    centroided(mse) <- FALSE  ## Takes quite some time.
    centroided(odmse) <- FALSE
    expect_identical(centroided(mse), centroided(odmse))
    ## Check error
    expect_that(centroided(odmse) <- c(TRUE, FALSE, TRUE), throws_error())

    ## peaksCount
    system.time(
        pk <- peaksCount(mse)
    )  ## 0.049
    system.time(
        pk2 <- peaksCount(odmse)
    )  ## 0.002
    expect_identical(pk, pk2)
    ## There's however something different if we're using removePeaks or clean.
    mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram", removePeaks=100000, clean=TRUE)
    odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk", removePeaks=100000, clean=TRUE)
    ## peaksCount
    system.time(
        pk <- peaksCount(mseRemPeaks)
    )  ## 0.046
    system.time(
        pk2 <- peaksCount(odmseRemPeaks)
    )  ## 51 secs.
    expect_identical(pk, pk2)

    ## spectra
    system.time(
        sp1 <- spectra(mse)
    )  ## 0.004
    sp1 <- lapply(sp1, function(z){
        z@polarity <- -1L
        return(z)
    })
    system.time(
        sp2 <- spectra(odmse)
    )  ## 29
    expect_identical(sp1, sp2)

    ## With clean and stuff.
    system.time(
        sp3 <- spectra(mseRemPeaks)
    )
    sp3 <- lapply(sp3, function(z){
        z@polarity <- -1L
        return(z)
    })
    system.time(
        sp4 <- spectra(odmseRemPeaks)
    )  ## 50.8
    expect_identical(sp3, sp4)
    ## do the queue afterwards.
    system.time(
        sp2 <- lapply(sp2, function(z){
            z <- MSnbase:::execute(odmseRemPeaks@spectraProcessingQueue[[1]], z)
            z <- MSnbase:::execute(odmseRemPeaks@spectraProcessingQueue[[2]], z)
            return(z)
        })
    ) ## 20 secs.
    expect_identical(sp2, sp4)

    ## Check some internal stuff....
    fd <- fData(odmse)
    fd <- fd[fd$fileIdx == 1, ]
    ## Get me the Spectrum1 objects.
    system.time(
        Test1 <- MSnbase:::.applyFun2SpectraOfFile(fData=fd, filenames=fileNames(odmse))
    ) ## 14.5 sec
    system.time(
        Test2 <- MSnbase:::.applyFun2SpectraOfFile2(fData=fd, filenames=fileNames(odmse))
    ) ## 3.7 sec
    expect_identical(Test1, Test2)
})


