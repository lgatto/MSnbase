readMSData <- function(files,
                       pdata = NULL,
                       msLevel = 2,
                       verbose = TRUE,
                       centroided = FALSE,
                       smoothed = FALSE,
                       removePeaks = 0,
                       clean = FALSE,
                       cache = 1) {
  ## TODO: add also a trimMz argument.
  if (msLevel == 1) ## cache currently only works for MS2 level data
    cache <- 0 
  msLevel <- as.integer(msLevel)
  if (!(msLevel > 0))
    stop("msLevel should be an integer > 0.")
  if (length(files) < 1)
    stop("At least one MS file is required.")
  extensions <- unique(toupper(sub("^.+\\.", "", files)))
  if (length(extensions) > 1)
    warning(paste("Reading different file formats in.",
                  "This is untested but you are welcome to try it out.",
                  "Please report back!", sep="\n"))
  ## Creating environment with Spectra objects
  assaydata <- new.env()
  for (f in files) {
    filen <- match(f, files) 
    msdata <- mzR::openMSfile(f)
    fullhd <- mzR::header(msdata)
    ## 
    if (msLevel == 1) {
      spidx <- which(fullhd$msLevel == 1)
      if (length(spidx)==0)
        stop("No MS1 spectra in file",f)
      if (verbose) {
        cat("Reading ", length(spidx), " MS1 spectra from file ",
            basename(f),"\n",sep="")
        pb <- txtProgressBar(min=0, max=length(spidx), style=3)
      }
      for (i in 1:length(spidx)) {
        if (verbose) setTxtProgressBar(pb, i)
        j <- spidx[i]
        hd <- fullhd[j,]
        sp <- new("Spectrum1",
                  peaksCount = hd$peaksCount,
                  rt = hd$retentionTime,
                  acquisitionNum = hd$acquisitionNum,
                  mz = mzR::peaks(msdata,j)[,1],
                  intensity = mzR::peaks(msdata,j)[,2],
                  fromFile = filen,
                  centroided = centroided)
        if (removePeaks > 0)
          sp <- removePeaks(sp, t=removePeaks)
        if (clean)
          sp <- clean(sp)
        assign(paste("X",i,sep=""), sp, assaydata)        
      }
    } else {
      spidx <- which(fullhd$msLevel > 1)
      if (length(spidx)==0)
        stop("No MS(n>1) spectra in file",f)
      if (verbose) {
        cat("Reading ",length(spidx)," MS2 spectra from file ",
            basename(f),"\n",sep="")
        pb <- txtProgressBar(min=0,max=length(spidx),style=3)
      }
      ## this was fullhd$acquisitionNum -- check/wrong 
      ## ms1scanNums <- MSnbase:::getBins(fullhd$msLevel[spidx])
      scanNums <- fullhd[fullhd$msLevel == 2,"precursorScanNum"]
      if (length(scanNums)!=length(spidx))
        stop("Number of spectra and precursor scan number do not match!")
      for (i in 1:length(spidx)) {
        if (verbose) setTxtProgressBar(pb, i)
        j <- spidx[i]
        hd <- fullhd[j,]
        sp <- new("Spectrum2",
                  precScanNum=as.integer(scanNums[i]),
                  precursorMz=hd$precursorMZ,
                  precursorIntensity=hd$precursorIntensity,
                  precursorCharge=hd$precursorCharge,
                  collisionEnergy=hd$collisionEnergy,
                  peaksCount=hd$peaksCount,
                  rt=hd$retentionTime,
                  acquisitionNum=hd$acquisitionNum,
                  mz=mzR::peaks(msdata,j)[,1],
                  intensity=mzR::peaks(msdata,j)[,2],
                  fromFile=filen,
                  centroided=centroided)
        if (removePeaks > 0)
          sp <- removePeaks(sp, t=removePeaks)
        if (clean)
          sp <- clean(sp)
        assign(paste("X",i,sep=""),sp,assaydata)
      }
    }
    if (verbose)
      close(pb)
    ## cache levels 2 and 3 not yet implemented
    cache <- testCacheArg(cache, maxCache=1)
    .cacheEnv <- setCacheEnv(assaydata, cache, lock=TRUE)
    ## if cache==2, do not lock,
    ## assign msdata in .cacheEnv
    ## then lock it
    ## and do not close(msdata) below; rm(msdata) is OK
    gc() ## could this help with Error in function (x): no function to return from, jumping to top level)...
    mzR::close(msdata)
    rm(msdata)
  }
  ## Create 'MSnProcess' object
  process <- new("MSnProcess",
                 processing=paste("Data loaded:",date()),
                 files=files,
                 smoothed=smoothed)  
  if (removePeaks > 0) {
    process@processing <- c(process@processing,
                            paste("Curves <= ", removePeaks, " set to '0': ", date(), sep=""))
  } else {
    if (clean)
      process@processing <- c(process@processing,
                              paste("Spectra cleaned: ",date(),sep=""))
  }
  ## Create 'fdata' and 'pdata' objects
  nms <- ls(assaydata)
  if (is.null(pdata)) {
    pdata <- new("NAnnotatedDataFrame",
                 dimLabels=c("sampleNames", "sampleColumns"))
  }
  fdata <- new("AnnotatedDataFrame",
               data=data.frame(
                 spectrum=1:length(nms),
                 row.names=nms))
  fdata <- fdata[ls(assaydata)] ## reorder features
  ## Create and return 'MSnPeaks' object
  if (verbose)
    cat("Creating 'MSnExp' object\n")
  toReturn <- new("MSnExp",
                  assayData = assaydata,
                  phenoData = pdata,
                  featureData = fdata,
                  processingData = process,
                  .cache = .cacheEnv)
  if (validObject(toReturn))
    return(toReturn)
}

