setMethod("writeMgfData",
          signature=signature("Spectrum"),
          function(object,
                   filename="spectrum.mgf",
                   COM=NULL,
                   TITLE=NULL) {
            if (file.exists(filename)) 
              message("Overwriting ",filename,".")
            if (is.null(COM))
              cat("COM=Spectrum exported by MSnbase on ",
                  date(),"\n", sep="", file=filename)
            writeMgfContent(object,TITLE=TITLE,filename=filename)
          })

setMethod("writeMgfData",
          signature=signature("MSnExp"),
          function(object,
                   filename="experiment.mgf",
                   COM=NULL) {
            if (file.exists(filename)) 
              warning("Overwriting ",filename,".")
            if (is.null(COM))
              cat("COM=Experiment exported by MSnbase on ",
                  date(), "\n", sep="", file=filename)
            tmp <- lapply(spectra(object),writeMgfContent,TITLE=NULL,filename)
          })

writeMgfContent <- function(sp,TITLE=NULL,filename) {
  cat("BEGIN IONS\n",file=filename,append=TRUE)  
  cat("SCANS=",acquisitionNum(sp),"\n",sep="",file=filename,append=TRUE)
  if (is.null(TITLE))
    cat("TITLE=MS",msLevel(sp),"spectrum\n",sep="",file=filename,append=TRUE)
  else 
    cat("TITLE=",TITLE,"\n",sep="",file=filename,append=TRUE)
  cat("RTINSECONDS=",rtime(sp),"\n",sep="",file=filename,append=TRUE)
  cat("PEPMASS=",precursorMz(sp),"\n",sep="",file=filename,append=TRUE)
  if (!is.na(precursorCharge(sp)))
    cat("CHARGE=",precursorCharge(sp),"+\n",sep="",file=filename,append=TRUE)
  dfr <- as(sp,"data.frame")
  tmp <- apply(dfr,1,
               function(x) cat(as.character(x),"\n",file=filename,append=TRUE))
  cat("END IONS\n",file=filename,append=TRUE)
}


## Code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
readMgfData <- function(file,
                        pdata = NULL,
                        centroided = TRUE,
                        smoothed = FALSE,
                        verbose = TRUE,
                        cache = 1) {
  if (verbose)
    cat("Scanning ",file,"...\n",sep="")
  mgf <- scan(file=file, what="", sep="\n", quote="", allowEscapes=FALSE, quiet=TRUE)
  begin <- grep("BEGIN IONS", mgf)
  end <- grep("END IONS", mgf)
  if (verbose) {
    cnt <- 1
    pb <- txtProgressBar(min = 0, max = length(begin), style = 3)
    
  }
  spectra <- vector("list",length=length(begin))
  fdata <- c()
  for (i in 1:length(begin)) {
    if (verbose) {
      setTxtProgressBar(pb, cnt)
      cnt <- cnt + 1
    }
    spec <- mgf[seq(begin[i], end[i])]
    spectra[[i]] <- mgfToSpectrum2(spec,centroided=centroided)
    ## adding the spectrum header as fData
    tmpdesc <- spec[grep("=",spec)]
    tmpdesc <- sapply(tmpdesc,strsplit,"=")
    desc <- unlist(lapply(tmpdesc,"[",2))
    names(desc) <- unlist(lapply(tmpdesc,"[",1))
    fdata <- rbind(fdata,desc)
  }
  if (verbose) 
    close(pb)
  nms <- paste("X",seq_len(length(spectra)),sep="")
  names(spectra) <- nms
  assaydata <- list2env(spectra)
  process <- new("MSnProcess",
                 processing=paste("Data loaded:",date()),
                 files=file,
                 smoothed=smoothed)
  if (is.null(pdata)) {
    pdata <- new("NAnnotatedDataFrame",
                 dimLabels=c("sampleNames", "sampleColumns"))
  }
  ## fdata <- new("AnnotatedDataFrame",
  ##              data=data.frame(spectrum=1:length(spectra),
  ##                row.names=nms))
  colnames(fdata) <- names(desc)
  rownames(fdata) <- nms
  fdata <- AnnotatedDataFrame(data=data.frame(fdata))
  fdata <- fdata[ls(assaydata)] ## reorder features
  ## only levels 0 and 1 for mgf peak lists
  cache <- testCacheArg(cache,maxCache=1)
  .cacheEnv <- newCacheEnv(assaydata, cache, lock=TRUE)
  toReturn <- new("MSnExp",
                  assayData = assaydata,
                  phenoData = pdata,
                  featureData = fdata,
                  processingData = process,
                  .cache = .cacheEnv)
  if (validObject(toReturn))
    return(toReturn)
}


mgfToSpectrum2 <- function(mgf, centroided) {
    mgf <- mgf[c(-1, -length(mgf))] ## remove "BEGIN IONS" and "END IONS"
    mgf <- sub("\t"," ",mgf) ## expecting mz and int to be separated by 1 space
    desc.idx <- grep("=",mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]
    ## mz.i <- adply(.data=mgf[-desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split=" ")))
    ## mz <- as.numeric(mz.i$V1)
    ## int <- as.numeric(mz.i$V2)
    mz.i <- sapply(spec,strsplit," ")
    mz <- as.numeric(unlist(lapply(mz.i,"[",1)))
    int <- as.numeric(unlist(lapply(mz.i,"[",2)))
    ## desc <- adply(.data=mgf[desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split="=")))
    ## desc <- desc[,-1]
    tmpdesc <- sapply(desc,strsplit,"=")
    desc <- unlist(lapply(tmpdesc,"[",2))
    names(desc) <- unlist(lapply(tmpdesc,"[",1))
    ## sp <- new("Spectrum2",
    ##           rt=as.numeric(desc[desc[,1] == "RTINSECONDS",2]),
    ##           acquisitionNum=as.integer(desc[desc[,1] == "SCANS",2]),
    ##           precursorMz=as.numeric(desc[desc[,1] == "PEPMASS",2]),
    ##           #precursorIntensity="",  ##I don't see any definition in MGF.
    ##           precursorCharge=as.integer(sub("[+,-]","",desc[desc[,1] == "CHARGE",2])),
    ##           #scanIndex="",
    ##           #collisionEnergy="",
    ##           peaksCount=length(int),
    ##           mz=mz,
    ##           intensity=int,
    ##           fromFile=1L,
    ##           centroided=centroided)
    sp <- new("Spectrum2",
              rt = ifelse("RTINSECONDS" %in% names(desc),as.numeric(desc["RTINSECONDS"]),0),
              ## acquisitionNum = ifelse("SCANS" %in% names(desc),as.integer(desc["SCANS"]),0L),
              precursorMz = ifelse("PEPMASS" %in% names(desc),as.numeric(desc["PEPMASS"]),0L),
              precursorCharge = ifelse("CHARGE" %in% names(desc),as.integer(sub("[+,-]","",desc["CHARGE"])),0L),
              peaksCount = length(int),
              mz = mz,
              intensity = int,
              fromFile = 1L,
              centroided = centroided)

    if(validObject(sp))
        return(sp)
}
