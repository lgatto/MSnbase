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
  cat("BEGIN IONS\n",file=filename,append=FALSE)  
  cat("SCANS=",acquisitionNum(sp),"\n",sep="",file=filename,append=TRUE)
  if (is.null(TITLE))
    cat("TITLE=MS ",msLevel(sp)," spectrum\n",sep="",file=filename,append=TRUE)
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
                        pdata=NULL,
                        centroided=TRUE,
                        smoothed=FALSE,
                        verbose=TRUE) {
    mgf <- scan(file=file, what="", sep="\n", quote="", allowEscapes=FALSE, quiet=TRUE)
    begin <- grep("BEGIN IONS", mgf)
    end <- grep("END IONS", mgf)
    if (verbose) {
      ._cnt <- 1
      pb <- txtProgressBar(min = 0, max = length(begin), style = 3)

    }
    spectra <- lapply(1:length(begin),
                      function(i) {
                        if (verbose) {
                          setTxtProgressBar(pb, ._cnt)
                          ._cnt <<- ._cnt+1
                        }
                        mgfToSpectrum2(mgf[seq(begin[i], end[i])],centroided=centroided)
                      })
    if (verbose) {
      close(pb)
      rm(pb)
      rm(._cnt)
    }
    nms <- paste("X",seq_len(length(spectra)),sep="")
    names(spectra) <- nms
    spectra.env <- list2env(spectra)
    process <- new("MSnProcess",
                   processing=paste("Data loaded:",date()),
                   files=file,
                   smoothed=smoothed)
    if (is.null(pdata)) {
        pdata <- new("NAnnotatedDataFrame",
                     dimLabels=c("sampleNames", "sampleColumns"))
    }
    fdata <- new("AnnotatedDataFrame",
                 data=data.frame(spectrum=1:length(spectra),
                 row.names=nms))
    fdata <- fdata[ls(spectra.env)] ## reorder features
    toReturn <- new("MSnExp",
                    assayData = spectra.env,
                    phenoData = pdata,
                    featureData = fdata,
                    processingData = process)
    if (validObject(toReturn))
        return(toReturn)
}


mgfToSpectrum2 <- function(mgf, centroided) {
    mgf <- mgf[c(-1, -length(mgf))] ## remove "BEGIN IONS" and "END IONS"
    desc.idx <- grep("=",mgf)
    mz.i <- adply(.data=mgf[-desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split=" ")))
    mz <- as.numeric(mz.i$V1)
    int <- as.numeric(mz.i$V2)
    desc <- adply(.data=mgf[desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split="=")))
    desc <- desc[,-1]
    sp <- new("Spectrum2",
              rt=as.numeric(desc[desc[,1] == "RTINSECONDS",2]),
              acquisitionNum=as.integer(desc[desc[,1] == "SCANS",2]),
              precursorMz=as.numeric(desc[desc[,1] == "PEPMASS",2]),
              #precursorIntensity="",  ##I don't see any definition in MGF.
              precursorCharge=as.integer(sub("[+,-]","",desc[desc[,1] == "CHARGE",2])),
              #scanIndex="",
              #collisionEnergy="",
              peaksCount=length(int),
              mz=mz,
              intensity=int,
              centroided=centroided)
    if(validObject(sp))
        return(sp)
}
