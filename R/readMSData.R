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
  if (all(unique(files) != files))
    stop("Non unique files provided as input. ")
  extensions <- unique(toupper(sub("^.+\\.", "", files)))
  if (length(extensions) > 1)
    warning(paste("Reading different file formats in.",
                  "This is untested and you are welcome to try it out.",
                  "Please report back!", sep="\n"))
  ## Creating environment with Spectra objects
  assaydata <- new.env()
  ioncount <- c()
  filenams <- filenums <- c()
  fullhd2 <- fullhdorder <- c()
  .instrumentInfo <- list()
  for (f in files) {
    filen <- match(f, files)
    filenums <- c(filenums, filen)
    filenams <- c(filenams, f)
    msdata <- mzR::openMSfile(f)    
    .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
    fullhd <- mzR::header(msdata)    
    ifelse(msLevel == 1,
           spidx <- which(fullhd$msLevel == 1),
           spidx <- which(fullhd$msLevel > 1))
    ## MS1 level
    if (msLevel == 1) {
      if (length(spidx) == 0)
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
                  tic = hd$totIonCurrent,
                  centroided = centroided)
        ioncount <- c(ioncount, sum(mzR::peaks(msdata,j)[,2]))
        if (removePeaks > 0)
          sp <- removePeaks(sp, t=removePeaks)
        if (clean)
          sp <- clean(sp)
        .fname <- paste0("X", i, ".", filen)
        assign(.fname, sp, assaydata)
        fullhdorder <- c(fullhdorder, .fname)
      }
    } else { ## MS>2 levels
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
                  tic = hd$totIonCurrent,                  
                  peaksCount=hd$peaksCount,
                  rt=hd$retentionTime,
                  acquisitionNum=hd$acquisitionNum,
                  mz=mzR::peaks(msdata,j)[,1],
                  intensity=mzR::peaks(msdata,j)[,2],
                  fromFile=filen,
                  centroided=centroided)
        ioncount <- c(ioncount, sum(mzR::peaks(msdata,j)[,2]))
        if (removePeaks > 0)
          sp <- removePeaks(sp, t=removePeaks)
        if (clean)
          sp <- clean(sp)
        .fname <- paste0("X", i, ".", filen)
        assign(.fname, sp, assaydata)
        fullhdorder <- c(fullhdorder, .fname)
      }
    }
    if (cache >= 1)
      fullhd2 <- rbind(fullhd2, fullhd[spidx, ]) 
    if (verbose)
      close(pb)
    gc() ## could this help with Error in function (x): no function to return from, jumping to top level)...
    mzR::close(msdata) ## DO NOT CLOSE IF CACHE LEVEL >= 2
    rm(msdata)
  }
  ## cache levels 2 and 3 not yet implemented
  cache <- testCacheArg(cache, maxCache = 1)
  if (cache >= 1) {
    fl <- sapply(assaydata, fromFile)
    featnms <- ls(assaydata) ## feautre names in final MSnExp
    fl <- fl[featnms] ## reorder file numbers
    stopifnot(all(sort(featnms) == sort(fullhdorder)))
    fullhdorder <- match(featnms, fullhdorder)    
    tmphd <- fullhd2[fullhdorder, ] ## reorder
    ioncount <- ioncount[fullhdorder]
    newhd <- data.frame(file = fl, 
                        retention.time = tmphd$retentionTime,
                        precursor.mz = tmphd$precursorMZ,
                        precursor.intensity = tmphd$precursorIntensity,
                        charge = tmphd$precursorCharge,
                        peaks.count = tmphd$peaksCount,
                        tic = tmphd$totIonCurrent,
                        ionCount = ioncount, 
                        ms.level = tmphd$msLevel,
                        acquisition.number = tmphd$acquisitionNum,
                        collision.energy = tmphd$collisionEnergy)
  } else {
    newhd <- NULL ## not used anyway
  }
  if (verbose)
    message("Caching...")
  .cacheEnv <- setCacheEnv(list("assaydata" = assaydata,
                                "hd" = newhd),
                           cache, lock = TRUE)
  ## if cache==2, do not lock,
  ## assign msdata in .cacheEnv
  ## then lock it
  ## and do not close(msdata) above; rm(msdata) is OK
  ## Create 'MSnProcess' object
  process <- new("MSnProcess",
                 processing = paste("Data loaded:",date()),
                 files = files,
                 smoothed = smoothed)  
  if (removePeaks > 0) {
    process@processing <- c(process@processing,
                            paste("Curves <= ", removePeaks, " set to '0': ", date(), sep=""))
  } else {
    if (clean)
      process@processing <- c(process@processing,
                              paste("Spectra cleaned: ", date(), sep=""))
  }
  ## Create 'fdata' and 'pdata' objects
  nms <- ls(assaydata)
  if (is.null(pdata)) {
    .pd <- data.frame(sampleNames = filenams,
                      fileNumbers = filenums)
    pdata <- new("NAnnotatedDataFrame",
                 data = .pd)
  }
  fdata <- new("AnnotatedDataFrame",
               data=data.frame(
                 spectrum=1:length(nms),
                 row.names=nms))
  fdata <- fdata[ls(assaydata)] ## reorder features
  ## expriment data slot
  if (length(.instrumentInfo) > 1) {
    cmp <- sapply(.instrumentInfo[-1], function(x) identical(x, .instrumentInfo[[1]]))
    if (!all(cmp)) {
      warning("According to the instrument information in the files, the data has been acquired on different instruments!")
      .instrumentInfo[[1]] <- list(manufacturer = paste(sapply(.instrumentInfo, "[[", "manufacturer"), collapse = ", "),
                                   model = paste(sapply(.instrumentInfo, "[[", "model"), collapse = ", "),
                                   ionisation = paste(sapply(.instrumentInfo, "[[", "ionisation"), collapse = ", "),
                                   analyzer = paste(sapply(.instrumentInfo, "[[", "analyzer"), collapse = ", "),
                                   detector = paste(sapply(.instrumentInfo, "[[", "detector"), collapse = ", "))
    } 
  } 
  expdata <- new("MIAPE",
                 instrumentManufacturer = .instrumentInfo[[1]]$manufacturer,
                 instrumentModel = .instrumentInfo[[1]]$model,
                 ionSource = .instrumentInfo[[1]]$ionisation,
                 analyser = .instrumentInfo[[1]]$analyzer,
                 detectorType = .instrumentInfo[[1]]$detector)
  ## Create and return 'MSnExp' object
  if (verbose)
    cat("Creating 'MSnExp' object\n")
  toReturn <- new("MSnExp",
                  assayData = assaydata,
                  phenoData = pdata,
                  featureData = fdata,
                  processingData = process,
                  experimentData = expdata,
                  .cache = .cacheEnv)
  if (validObject(toReturn))
    return(toReturn)
}

