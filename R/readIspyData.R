##' Reads an ispy2 result spread sheet and creates a fully featured MSnSet instance
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

