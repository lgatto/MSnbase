context("OnDiskMSnExp class, other methods")

############################################################
## Load the required data files.
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
## Load the data as an MSnExp into memory.
mse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="ram")
## Load the data as OnDiskMSnExp.
odmse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="disk")
## All the same with removePeaks.
mseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="ram",
                          removePeaks=10000, clean=TRUE)
odmseRemPeaks <- readMSData(files=mzf, msLevel=1, backend="disk",
                            removePeaks=10000, clean=TRUE)


############################################################
## plot
test_that("OnDiskMSnExp plot", {
    ## Would be nice to know what the plot function is actually doing though...
    ## seems I can forget that for larger experiments; takes way to long.
})

############################################################
## trimMz
test_that("OnDiskMSnExp trimMz", {
    ## Comparing timings and results for the trimMz.
    system.time(
        mseT <- trimMz(mse, mzlim=c(300, 310))
    ) ## 7.3 sec.
    system.time(
        odmseT <- trimMz(odmse, mzlim=c(300, 310))
    ) ## woah, 0.009 sec (what a surprise ;) )
    ## Test the results.
    system.time(
        mseTmz <- mz(mseT)
    ) ## 0.04 sec
    system.time(
        odmseTmz <- mz(odmseT)
    ) ##  sec
    expect_identical(mseTmz, odmseTmz)
    system.time(
        mseTsp <- spectra(mseT)
    ) ## 0.005
    system.time(
        odmseTsp <- spectra(odmseT)
    ) ## 6.8
    odmseTsp <- lapply(odmseTsp, function(z){
        z@polarity <- integer()
        z@scanIndex <- integer()
    })
    expect_identical(mseTmz, odmseTmz)
})

