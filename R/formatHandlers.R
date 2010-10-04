##' Load mzXML data as an MSnExp object
##'
##' This function reads the MS1 or MS2 spectra defined in an
##' mzXML file and stores them, as well as miscellaneous
##' meta-data, in an object of class \code{"\linkS4class{MSnExp}"}.
##' Rudimentary data processing like low intensity peaks and
##' zero intensity data points removal can optionally already be
##' applied at this stage.
##' 
##' @title Load mzXML data
##' @usage readMzXMLData(files,msLevel,verbose,centroided,smoothed,description,removePeaks,clean)
##' @param files a character of mzXML file(s) to be read.
##' @param msLevel a numeric indicating the MS level of spectra to be
##' imported (currently 1 or 2 only), default is 2.
##' @param verbose a logical value indication whether a progress should
##' report progress of data acquisition.
##' @param centroided a logical value indicating if peaks are centroided,
##' default is FALSE.
##' @param smoothed a logical value indicating if peak curves are smoothed,
##' default is FALSE.
##' @param description a character string describing the experiment. 
##' @param removePeaks a numeric value to define low-intensity peaks
##' that should be removed; indicate 0 (default) to keep data as is.
##' @param clean a logical value specifying whether MZ data points
##' with 0 intensity (if for instance removePeaks > 1) should be
##' completely removed.
##' @return an object of class \code{"\linkS4class{MSnExp}"} 
##' @seealso see \code{"\linkS4class{MSnExp}"} class as well as
##' removePeaks and clean methods for rudimentary data processing
##' after being loaded.
##' @TODO Add links to removePeaks and clean methods once documented.
##' @examples
##' mzxmlfile <- dir(system.file("inst/extdata",package="MSnbase"),pattern="mzXML$")
##' msnexp <- readMzXMLData(mzxmlfile)
##' msnexp
##' @author Laurent Gatto
##' @export 
##' @keywords data utilities 
readMzXMLData <- function(files,
                          msLevel=2,
                          verbose=TRUE,
                          centroided=FALSE,
                          smoothed=FALSE,
                          description=character(),
                          removePeaks=0,
                          clean=FALSE) {
  msLevel <- match.arg(msLevel,1:2,several.ok=FALSE)
  ## Opening file handles to mzXML files
  n <- length(files)
  if (n<1)
    stop("At least one mzXML file is required.")
  id <- numeric()
  for (i in 1:n) {
    if (!xcms:::rampIsFile(files[i]))
      stop("No such mzXML file: ",files[i])
    id[i] <- xcms:::rampOpen(files[i])
  }
  on.exit(xcms:::rampCloseAll())
  on.exit(gc(),add=TRUE)
  
  ## Acquiring peak data and create 'MSnPeaksData' object  
  fromFile <- c()
  spectra <- vector("list",length=n)
  for (i in 1:n) {
    if (msLevel>1) {
      raw <- xcms:::rampRawDataMSn(id[i])
    } else {
      raw <- xcms:::rampRawData(id[i])
    }
    if (verbose)
      cat("Acquiring data for",length(raw$rt),
          "precursors from file",files[i],"\n")
    spectra[[i]] <- vector("list",length=length(raw$rt))
    fromFile <- c(fromFile,rep(i,length(raw$rt)))
    pb <- txtProgressBar(min=0,max=length(raw$rt),style=3)
    for (j in 1:length(raw$rt)) {
      if (msLevel>1) {
        spectra[[i]][[j]] <- rawToSpectrum2(raw,j,
                                            clean=clean,
                                            removePeaks=removePeaks,
                                            updatePeaksCount=TRUE,
                                            verbose=FALSE)
      } else {
        spectra[[i]][[j]] <- rawToSpectrum1(raw,j,
                                            clean=clean,
                                            removePeaks=removePeaks,
                                            verbose=FALSE)
      }
      if (verbose)
        setTxtProgressBar(pb, j)
    }
  close(pb)
  }
  rm(raw)
  gc()
  
  ## Create 'MSnPeaksPocessing' object
  process <- new("MSnProcess",
                 processing=paste("Data Loaded:",date()),
                 smoothed=smoothed,
                 centroided = centroided)

  if (removePeaks>0) {
    process@processing <- c(process@processing,
                            paste("Curves <= ",removePeaks," set to '0': ",date(),sep=""),
                            paste("Spectra cleaned: ",date(),sep=""))
  } else {
    if (clean)
      process@processing <- c(process@processing,
                              paste("Spectra cleaned: ",date(),sep=""))
  }
  spectra <- unlist(spectra)
  if (msLevel>1) {
    ms1scanNums <- getBins(sapply(spectra,acquisitionNum))
    if (length(ms1scanNums)!=length(spectra))
      stop("Number of spectra and ms1scan numbers do not match!")
    for (i in 1:length(spectra)) {
      spectra[[i]]@ms1scan <- as.integer(ms1scanNums[i])
    }
  }

  ## Create and return 'MSnPeaks' object
  if (verbose)
    cat("Creating 'MSnExp' object\n")
  toReturn <- new("MSnExp",
                  spectra=spectra,
                  process=process,
                  description=description,
                  fromFile=fromFile,
                  files=files)
  return(toReturn)
}



## readMzXMLData.old <- function(files,
##                               verbose=TRUE,
##                               centroided=FALSE,
##                               smoothed=FALSE,
##                               description=character(),
##                               removePeaks=0,
##                               clean=FALSE) {
##   ## Opening file handles to mzXML files
##   n <- length(files)
##   id <- numeric()
##   for (i in 1:n) {
##     if (!xcms:::rampIsFile(files[i]))
##       stop("No such mzXML file: ",files[i])
##     id[i] <- xcms:::rampOpen(files[i])
##   }
##   on.exit(xcms:::rampCloseAll())
##   on.exit(gc(),add=TRUE)  
##   ## Acquiring peak data and create 'MSnPeaksData' object  
##   fromFile <- c()
##   spectra <- vector("list",length=n)  
##   for (i in 1:n) {
##     raw <- xcms:::rampRawDataMSn(id[i])
##     peaks <- list(msnRt = raw$rt, 
##                   msnAcquisitionNum = raw$acquisitionNum,
##                   msnPrecursorNum = raw$precursorNum,
##                   msnPrecursorMz = raw$precursorMZ,
##                   msnPrecursorIntensity = raw$precursorIntensity,
##                   msnPeaksCount = raw$peaksCount,
##                   msnLevel = raw$msLevel,
##                   msnPrecursorCharge = raw$precursorCharge,
##                   msnScanindex = raw$scanindex,
##                   msnCollisionEnergy = raw$collisionEnergy,
##                   msnMz = raw$mz,
##                   msnIntensity = raw$intensity)
##     if (verbose)
##       cat("Acquiring data for",length(raw$rt),
##           "precursors from file",files[i],"\n")
##     spectra[[i]] <- vector("list",length=length(raw$rt))
##     fromFile <- c(fromFile,rep(i,length(raw$rt)))
##     pb <- txtProgressBar(min=0,max=length(raw$rt),style=3)
##     for (j in 1:length(raw$rt)) {
##       spectra[[i]][[j]] <- rawToSpectrum2(peaks,j,
##                                           clean=clean,
##                                           removePeaks=removePeaks,
##                                           verbose=FALSE)
##       if (verbose)
##         setTxtProgressBar(pb, j)
##     }
##   close(pb)
##   }
##   rm(raw)
##   gc()  
##   ## Create 'MSnPeaksPocessing' object
##   process <- new("MSnProcess",
##                  processing=paste("Data Loaded:",date()),
##                  smoothed=smoothed,
##                  centroided = centroided)
##   if (removePeaks>0) {
##     process@processing <- c(process@processing,
##                             paste("Curves <= ",removePeaks," set to '0': ",date(),sep=""),
##                             paste("Spectra cleaned: ",date(),sep=""))
##   } else {
##     if (clean)
##       process@processing <- c(process@processing,
##                               paste("Spectra cleaned: ",date(),sep=""))
##   }
##   spectra <- unlist(spectra)
##   ms1scanNums <- getBins(sapply(spectra,acquisitionNum))
##   if (length(ms1scanNums)!=length(spectra))
##     stop("Number of spectra and ms1scan numbers do not match!")
##   for (i in 1:length(spectra)) {
##     spectra[[i]]@ms1scan <- as.integer(ms1scanNums[i])
##     }
##   ## Create and return 'MSnPeaks' object
##   if (verbose)
##     cat("Creating 'MSnExp' object\n")
##   toReturn <- new("MSnExp",
##                   spectra=spectra,
##                   process=process,
##                   description=description,
##                   fromFile=fromFile,
##                   files=files)
##   return(toReturn)
## }

