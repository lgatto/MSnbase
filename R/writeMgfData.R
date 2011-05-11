setMethod("writeMgfData",
          signature=signature("Spectrum"),
          function(object,
                   filename="spectrum.mgf",
                   COM=NULL,
                   TITLE=NULL) {
            if (file.exists(filename)) 
              message("Overwriting ",filename,".")
            if (is.null(COM))
              cat("COM=Spectrum exported by MSnbase on ",date(),"\n",sep="",file=filename)
            writeMgfContent(object,TITLE=TITLE,filename=filename)
          })

setMethod("writeMgfData",
          signature=signature("MSnExp"),
          function(object,
                   filename="experiment.mgf",
                   COM=NULL) {
            if (file.exists(filename)) 
              message("Overwriting ",filename,".")
            if (is.null(COM))
              cat("COM=Experiment exported by MSnbase on ",date(),"\n",sep="",file=filename)
            tmp <- lapply(spectra(object),writeMgfContent,TITLE=NULL,filename)
          })

writeMgfContent <- function(sp,TITLE=NULL,filename) {
  cat("BEGIN IONS\n",file=filename,append=TRUE)  
  cat("SCANS=",acquisitionNum(sp),"\n",sep="",file=filename,append=TRUE)
  if (is.null(TITLE))
    cat("TITLE=MS ",msLevel(sp)," sp\n",sep="",file=filename,append=TRUE)
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

