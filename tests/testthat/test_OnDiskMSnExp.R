context("OnDiskMSnExp class")

############################################################
## Testing the on-disk MSnExp stuff.
test_that("OnDiskMSnExp constructor", {
})

test_that("read and validate OnDiskMSnExp data", {
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
    ## Reading the stuff in memory.
    mse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="ram")
    expect_identical(as.character(class(mse)), "MSnExp")
    ## Read as an OnDiskMSnExp
    odmse <- readMSData(files=mzf, msLevel=1, centroided=TRUE, backend="disk")

    ## Check if we get the same data! MSnExp will retrieve that from the spectra,
    ## OnDiskMSnExp from featureData.
    ## fromFile
    expect_identical(fromFile(mse), fromFile(odmse))

    ## msLevel
    expect_identical(msLevel(mse), msLevel(odmse))
})


