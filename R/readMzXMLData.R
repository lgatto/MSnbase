readMzXMLData <- function(files,
                          pdata=NULL,
                          msLevel=2,
                          verbose=TRUE,
                          centroided=FALSE,
                          smoothed=FALSE,
                          removePeaks=0,
                          clean=FALSE) {
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
          "spectra from file",files[i],"\n")
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
      spectra[[i]][[j]]@fromFile <- i
      if (verbose)
        setTxtProgressBar(pb, j)
    }
    close(pb)
  }
  rm(raw)
  gc()
  
  ## Create 'MSnPocessing' object
  process <- new("MSnProcess",
                 processing=paste("Data loaded:",date()),
                 files=files,
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

  
  nms <- paste("X",seq_len(length(spectra)),sep="")
  names(spectra) <- nms

  spectra.env <- list2env(spectra)

  ## here, a multiplex number (possibly 1) should be used
  ## to generate th NAnnotatedDataFrame object to afterwards,
  ## when the corresponding MSnSet is instanciated,
  ## automatically create the AnnotatedDataFrame.
  if (is.null(pdata)) {
    pdata <- new("NAnnotatedDataFrame",
                 dimLabels=c("sampleNames", "sampleColumns"))
  }
  fdata <- new("AnnotatedDataFrame",
               data=data.frame(spectum=1:length(spectra),
                 row.names=nms))
  fdata <- fdata[ls(spectra.env)] ## reorder features
  
  ## Create and return 'MSnPeaks' object
  if (verbose)
    cat("Creating 'MSnExp' object\n")
  
  toReturn <- new("MSnExp",
                  assayData=spectra.env,
                  phenoData=pdata,
                  featureData=fdata,
                  process=process)
  return(toReturn)
}



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
##' @usage readMzXMLData(files,msLevel,verbose,centroided,smoothed,removePeaks,clean)
##' @param files a character of mzXML file(s) to be read.
##' @param msLevel a numeric indicating the MS level of spectra to be
##' imported (currently 1 or 2 only), default is 2.
##' @param verbose a logical value indication whether a progress should
##' report progress of data acquisition.
##' @param centroided a logical value indicating if peaks are centroided,
##' default is FALSE.
##' @param smoothed a logical value indicating if peak curves are smoothed,
##' default is FALSE.
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
##' @references
##' Information about the mzXML format as well converters from vendor
##' specific formats to mzXML: 
##' \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML}.
##' @note
##' A new open cummunity-developed format \code{mzML} (Martens
##' \texit{et al.}, published in MCP in 2010 --
##' \url{http://www.mcponline.org/content/early/2010/08/17/mcp.R110.000133.abstract}), 
##' developed under the banner of the HUPO PSI is expected
##' to supercede \code{mzXML} and will be implemented in the MSnbase.
##' @author Laurent Gatto
##' @export 
##' @keywords data utilities 
readMzXMLData.old <- function(files,
                              msLevel=2,
                              verbose=TRUE,
                              centroided=FALSE,
                              smoothed=FALSE,
                              removePeaks=0,
                              clean=FALSE) {
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
          "spectra from file",files[i],"\n")
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
                 processing=paste("Data loaded:",date()),
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
                  fromFile=fromFile,
                  files=files)
  return(toReturn)
}

