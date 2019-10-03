
.applyFun2IndividualSpectraOfFile <- function(fData, filenames,
                                                  queue = NULL,
                                                  APPLYFUN = NULL,
                                                  fastLoad = TRUE,
                                                  ...) {
    suppressPackageStartupMessages(
        require(MSnbase, quietly = TRUE)
    )
    verbose. <- isMSnbaseVerbose()
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    if (length(queue) > 0) {
        if (verbose.) {
            message("Apply lazy processing step(s):")
            for (j in 1:length(queue))
                message(" o '", queue[[j]]@FUN, "' with ",
                        length(queue[[j]]@ARGS), " argument(s).")
        }
    }
    filename <- filenames[fData$fileIdx[1]]
    ## issue #214: define backend based on file format.
    fileh <- MSnbase:::.openMSfile(filename)
    ## Reading the header for the selecte spectra. This is to avoid getting
    ## "memory not mapped" errors when reading mz and intensity values from
    ## certain mzML files (issue #170). Since this problem seems to be absent
    ## on linux and Windows systems we allow the user to disable it.
    ## Also we are just reading the header for the last spectrum since that
    ## seems to fix it too.
    if (!fastLoad)
        hd_spectra <- mzR::header(fileh, max(fData$spIdx))
    n_rows <- nrow(fData)
    do_queue <- length(queue)
    do_apply <- !is.null(APPLYFUN)
    res <- vector("list", n_rows)
    for (i in 1:n_rows) {
        pks <- mzR::peaks(fileh, fData$spIdx[i])
        if (fData$msLevel[i] == 1) {
            sp <- MSnbase:::Spectrum1_mz_sorted(rt = fData$retentionTime[i],
                                      acquisitionNum = fData$acquisitionNum[i],
                                      scanIndex = fData$spIdx[i],
                                      tic = fData$totIonCurrent[i],
                                      mz = pks[, 1],
                                      intensity = pks[, 2],
                                      fromFile = fData$fileIdx[i],
                                      centroided = fData$centroided[i],
                                      smoothed = fData$smoothed[i],
                                      polarity = fData$polarity[i])
        } else {
            sp <- MSnbase:::Spectrum2_mz_sorted(msLevel = fData$msLevel[i],
                                      rt = fData$retentionTime[i],
                                      acquisitionNum = fData$acquisitionNum[i],
                                      scanIndex = fData$spIdx[i],
                                      tic = fData$totIonCurrent[i],
                                      mz = pks[, 1],
                                      intensity = pks[, 2],
                                      fromFile = fData$fileIdx[i],
                                      centroided = fData$centroided[i],
                                      smoothed = fData$smoothed[i],
                                      polarity = fData$polarity[i],
                                      merged = fData$mergedScan[i],
                                      precScanNum = fData$precursorScanNum[i],
                                      precursorMz = fData$precursorMZ[i],
                                      precursorIntensity = fData$precursorIntensity[i],
                                      precursorCharge = fData$precursorCharge[i],
                                      collisionEnergy = fData$collisionEnergy[i])
        }
        ## And now go through the processing queue - if not empty...
        if (do_queue) {
            for (pStep in queue) {
                sp <- MSnbase::executeProcessingStep(pStep, sp)
            }
        }
        ## Apply the function, if provided
        if (do_apply)
            res[[i]] <- APPLYFUN(sp, ...)
        else
            res[[i]] <- sp
    }
    names(res) <- rownames(fData)
    mzR::close(fileh)
    rm(fileh)
    ## Intermediate #151 fix. Performance-wise would be nice to get rid of this.
    ## gc()
    res
}

## To use this along with profvis
## library(profvis)
## source("/Users/jo/Projects/git/MSnbase/inst/scripts/spectrapply_test.R")
## library(mzR)
## fl <- "/Users/jo/data/2017/2017_02/090217_132h_RT_a.mzML"
fh <- mzR::openMSfile(fl)

hdr <- mzR::header(fh)
mzR::close(fh)
fData <- hdr
fData$spIdx <- hdr$seqNum
fData$fileIdx <- 1L
fData$smoothed <- FALSE
fData$centroided <- TRUE

fastLoad <- FALSE


tmp <- profvis(.applyFun2IndividualSpectraOfFile(fData, fl,
                                                 fastLoad = fastLoad),
               interval = 0.005)

tmp <- profvis(.applyFun2IndividualSpectraOfFile(fData, fl,
                                                 fastLoad = fastLoad,
                                                 APPLYFUN = mz),
               interval = 0.005)
