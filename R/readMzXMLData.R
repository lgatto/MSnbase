readMzXMLData <- function(files,
                          pdata = NULL,
                          msLevel = 2,
                          verbose = TRUE,
                          centroided = FALSE,
                          smoothed = FALSE,
                          removePeaks = 0,
                          clean = FALSE,
                          cache = 1) {

  .Defunct(new="readMSData")
  ## msg <- c("The 'readMzXMLData' function is deprecated\n",
  ##         "Please use 'readMSData' instead.")
  ## .Deprecated(msg=msg)
  ## msLevel <- as.integer(msLevel)
  ## ## Opening file handles to mzXML files
  ## n <- length(files)
  ## if (n<1)
  ##   stop("At least one mzXML file is required.")
  ## id <- numeric()
  ## for (i in 1:n) {
  ##   if (!xcms:::rampIsFile(files[i]))
  ##     stop("No such mzXML file: ",files[i])
  ##   id[i] <- xcms:::rampOpen(files[i])
  ## }
  ## on.exit(xcms:::rampCloseAll())
  ## on.exit(gc(),add=TRUE)
  
  ## ## Acquiring peak data and create 'MSnPeaksData' object  
  ## spectra <- vector("list",length=n)
  ## for (i in 1:n) {
  ##   if (msLevel>1) {
  ##     raw <- xcms:::rampRawDataMSn(id[i])
  ##   } else {
  ##     raw <- xcms:::rampRawData(id[i])
  ##   }
  ##   if (verbose)
  ##     cat("Acquiring data for",length(raw$rt),
  ##         "spectra from file",files[i],"\n")
  ##   spectra[[i]] <- vector("list",length=length(raw$rt))
  ##   if (verbose)
  ##     pb <- txtProgressBar(min=0,max=length(raw$rt),style=3)
  ##   for (j in 1:length(raw$rt)) {
  ##     if (msLevel>1) {
  ##       spectra[[i]][[j]] <- rawToSpectrum2(raw,j,
  ##                                           clean=clean,
  ##                                           centroided=centroided,
  ##                                           removePeaks=removePeaks,
  ##                                           updatePeaksCount=TRUE,
  ##                                           verbose=FALSE)
  ##     } else {
  ##       spectra[[i]][[j]] <- rawToSpectrum1(raw,j,
  ##                                           clean=clean,
  ##                                           removePeaks=removePeaks,
  ##                                           centroided=centroided,
  ##                                           verbose=FALSE)
  ##     }
  ##     spectra[[i]][[j]]@fromFile <- i
  ##     if (verbose)
  ##       setTxtProgressBar(pb, j)
  ##   }
  ##   if (verbose)
  ##     close(pb)
  ## }
  ## rm(raw)
  ## gc()
  
  ## ## Create 'MSnPocessing' object
  ## process <- new("MSnProcess",
  ##                processing=paste("Data loaded:",date()),
  ##                files=files,
  ##                smoothed=smoothed)
  
  ## if (removePeaks>0) {
  ##   process@processing <- c(process@processing,
  ##                           paste("Curves <= ",removePeaks," set to '0': ",date(),sep=""),
  ##                           paste("Spectra cleaned: ",date(),sep=""))
  ## } else {
  ##   if (clean)
  ##     process@processing <- c(process@processing,
  ##                             paste("Spectra cleaned: ",date(),sep=""))
  ## }
  ## spectra <- unlist(spectra)
  ## if (msLevel>1) {
  ##   scanNums <- getBins(sapply(spectra,acquisitionNum))
  ##   if (length(scanNums)!=length(spectra))
  ##     stop("Number of spectra and precursor scan number do not match!")
  ##   for (i in 1:length(spectra)) {
  ##     spectra[[i]]@precScanNum <- 0L ## as.integer(ms1scanNums[i])
  ##   }
  ## }

  
  ## nms <- paste("X",seq_len(length(spectra)),sep="")
  ## names(spectra) <- nms

  ## assaydata <- list2env(spectra)

  ## ## here, a multiplex number (possibly 1) should be used
  ## ## to generate th NAnnotatedDataFrame object to afterwards,
  ## ## when the corresponding MSnSet is instanciated,
  ## ## automatically create the AnnotatedDataFrame.
  ## if (is.null(pdata)) {
  ##   pdata <- new("NAnnotatedDataFrame",
  ##                dimLabels=c("sampleNames", "sampleColumns"))
  ## }
  ## fdata <- new("AnnotatedDataFrame",
  ##              data=data.frame(spectrum=1:length(spectra),
  ##                row.names=nms))
  ## fdata <- fdata[ls(assaydata)] ## reorder features
  ## ## cache levels 2 and 3 not yet implemented
  ## cache <- testCacheArg(cache,maxCache=1)
  ## .cacheEnv <- setCacheEnv(assaydata, cache, lock=TRUE)
  ## ## Create and return 'MSnPeaks' object
  ## if (verbose)
  ##   cat("Creating 'MSnExp' object\n")
  ## toReturn <- new("MSnExp",
  ##                 assayData = assaydata,
  ##                 phenoData = pdata,
  ##                 featureData = fdata,
  ##                 processingData = process,
  ##                 .cache = .cacheEnv)
  ## if (validObject(toReturn))
  ##   return(toReturn)
}
