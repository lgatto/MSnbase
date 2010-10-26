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
readMzXMLData <- function(files,
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




##' Reads an ispy2 result spread sheet and creates a fully features MSnSet instance
##'
##'\code{readIspyData} reads an ispy2 result spread sheet in tsv format
##' and parses the reporter values and the identification meta data to
##' generate a fully features MSnSet object.
##' 
##' @title readIspyData
##' @aliases readIspyData
##' @usage
##' readIspyData(file,uniquePeps=TRUE,pep=0.05,
##' na.rm=TRUE, min.int=0,reporters=19:23,skipFullUp=24,verbose=TRUE) 
##' @param file A character indicating the name of ispy2 result tsv file
##' @param uniquePeps A logical indicating whether only unique peptides
##' should be imported. Default is TRUE.
##' @param pep A numeric indicating the posterior error probability thershold
##' for peptides to be considered correctly identified. Default is 0.05.
##' @param na.rm A logical indicating whether reporter ions containing one or
##' more NA values should be excluded. Default is TRUE.
##' @param min.int A numeric indicating the minimal summed intensity threshold
##' for reporter data to be imported. Default is 0.
##' @param reporters A numeric vector indicating the reporter columns.
##' Default is 19 to 23 (for 4 reporter ions).
##' @param skipFillUp A numeric indicating columns that should not be filled
##' up. Generally columns after the reporter ions (see reporters parameter).
##' Default is 24.
##' @param verbose A logical indicating whether verbose output is to be
##' printed out.
##' @return an object of class \code{\linkS4class{MSnSet}
##' @references 
##' @note ispy is a set of perl script to analyse MS1 and MSMS data
##' developed by Phil D. Charles <pdc35@cam.ac.uk> at CCP
##' \url{http://www.bio.cam.ac.uk/proteomics/}
##' @seealso \code{\linkS4class{MSnSet}} class.
##' @author Laurent Gatto
##' @export
##' @keywords data utilities 
readIspyData <- function(file="ispy_results.tsv",
                         uniquePeps=TRUE,
                         pep=0.05,
                         na.rm=TRUE,
                         min.int=0,
                         reporters=19:23,
                         skipFillUp=24,
                         verbose=TRUE) {
  .fillUp <- function(x) {
    ## Fills up the left-most columns of
    ## the ispy results spread sheet
    if (!any(is.na(x)) & !any(x!=""))
      return(x)
    for (i in 2:length(x)) {
      if (is.na(x[i])) x[i] <- x[i-1]
      if (x[i]=="") x[i] <- x[i-1]
    }
    return(x)
  }
  if (verbose)
    cat("Reading table\n")
  tab <- read.csv(file,header=TRUE,sep="\t",fill=TRUE,
                  colClasses=c(
                    "numeric",       # Index
                    "factor",        # ProteinAccession
                    "character",     # ProteinDescription
                    "numeric",       # ProteinQ.Value.allpeptides
                    "numeric",       # ProteinType1Error.allpeptides
                    "numeric",       # ProteinQ.Value.uniquepeptidesonly.
                    "numeric",       # ProteinType1Error.uniquepeptidesonly.
                    "numeric",       # NumberOfUniquePeptides.FDR.1pc
                    "character",     # PeptideSequence
                    "numeric",       # PeptideParentProteins
                    "character",     # MSFileID
                    "numeric",       # PeptideType1ErrorForMSFile
                    "factor",        # VariableModifications
                    "numeric",       # VariableModificationMassShift
                    "numeric",       # PrecursorMZ
                    "character",     # MSScanID
                    "character",     # MascotQuery
                    "numeric",       # PosteriorErrorProbability
                    rep("numeric",length(reporters)),
                    "numeric"))      # precursorrelativesignal
  names(tab) <- gsub("\\.\\.","pc",gsub("_","",names(tab)))
  names(tab)[ncol(tab)] <- "PrecursorRelativeSignal"
  tab$MascotQuery <- sub(",[0-0]+\\)&","",
                         sub("=HYPERLINK\\(","",tab$MascotQuery))
  .exprs <- as.matrix(tab[,reporters])
  .featureData <- tab[,c(-reporters,-skipFillUp)]
  if (verbose)
    cat("Filling up table\n")
  pb <- txtProgressBar(min = 0, max = ncol(.featureData), style = 3)
  for (i in 1:ncol(.featureData)) {
    .featureData[,i] <- .fillUp(.featureData[,i])
    setTxtProgressBar(pb,i)
  }
  close(pb)   
  if (verbose)
    cat("Filtering\n")
  ## Default filters keep all features
  keep.uniq <- keep.pep <- keep.na <- keep.int <- rep(TRUE,nrow(.exprs))
  ## Filtering on uniqueness of peptides
  if (uniquePeps) 
    keep.uniq <- .featureData$PeptideParentProteins==1
  ## Filtering on posterior error probability
  keep.pep <- .featureData$PosteriorErrorProbability<pep
  ## Remove features with NAs
  if (na.rm)
    keep.na <- apply(.exprs,1,function(x) !any(is.na(x)))
  rowsums <- rowSums(.exprs,na.rm=TRUE)
  keep.int <- rowsums>min.int
  if (verbose)
    cat(" keep.na: ",sum(keep.na),"\n",
        "keep.int: ",sum(keep.int),"\n",
        "keep.pep: ",sum(keep.pep),"\n",
        "keep.uniq: ",sum(keep.uniq),"\n")
  ## Applying filter
  keep <- keep.na & keep.pep & keep.uniq & keep.int
  .featureData <- .featureData[keep,]
  .exprs <- .exprs[keep,]
  if (any(is.na(.featureData)))
    stop("NA values in .featureData.")
  if (any(is.na(.exprs)))
    warning("NA values in .exprs.")
  ## Preparing return object slots and return new MSnSet
  .process <- new("MSnProcess",
                  processing=paste("Data loaded:",date(),
                    "from ispy result file."),
                  normalised=FALSE)
  return(new("MSnSet",
             exprs=.exprs,
             featureData=new("AnnotatedDataFrame",data=.featureData),
             process=.process,
             files=file))
}

