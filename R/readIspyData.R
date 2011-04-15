readIspyData <- function(file="ispy_results.tsv",
                         uniquePeps=TRUE,
                         pep=0.05,
                         na.rm=TRUE,
                         min.int=0,
                         reporters=19:23,
                         skipFillUp=24,
                         fillUp=FALSE,
                         keepAll=FALSE,
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
                    "numeric",       # ProteinQ.Value.uniquepeptidesonly.
                    "numeric",       # ProteinType1Error.uniquepeptidesonly.
                    "numeric",       # NumberOfUniquePeptides.FDR.1pc
                    "character",     # PeptideSequence
                    "numeric",       # PeptideParentProteins
                    "character",     # MSFileID
                    "numeric",       # PeptideType1ErrorForMSFile
                    "factor",        # Fixed_Modifications
                    "numeric",       # FixedModificationMassShift                    
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
  ## tab$MascotQuery <- sub(",[0-0]+\\)&","",
  ##                        sub("=HYPERLINK\\(","",tab$MascotQuery))
  tab$MascotQuery <- sub("^.+F0","F0",
                         sub("\\.dat.+$",".dat",tab$MascotQuery))
  .exprs <- as.matrix(tab[,reporters])
  .featureData <- tab[,-reporters] ## was c(-reporters,-skipFillUp)
  if (fillUp) {
    if (verbose)
      cat("Filling up table\n")
    if (verbose)
      pb <- txtProgressBar(min = 0, max = ncol(.featureData), style = 3)
    for (i in 1:ncol(.featureData)) {
      .featureData[,i] <- .fillUp(.featureData[,i])
      if (verbose)
        setTxtProgressBar(pb,i)
    }
    if (verbose)
      close(pb)
  }
  if (verbose)
    cat("Filtering\n")
  ## Default filters keep all features
  keep.uniq <- keep.pep <- keep.na <- keep.int <- rep(TRUE,nrow(.exprs))
  ## Filtering on uniqueness of peptides
  if (uniquePeps) 
    keep.uniq <- .featureData$PeptideParentProteins==1
  ## Filtering on posterior error probability
  keep.pep <- .featureData$PosteriorErrorProbability<=pep
  ## Remove features with NAs
  if (na.rm)
    keep.na <- apply(.exprs,1,function(x) !any(is.na(x)))
  rowsums <- rowSums(.exprs,na.rm=TRUE)
  keep.int <- rowsums>=min.int
  keep.int[is.na(keep.int)] <- TRUE
  if (verbose & !keepAll)
    cat(" keep.na: ",sum(keep.na),"\n",
        "keep.int: ",sum(keep.int),"\n",
        "keep.pep: ",sum(keep.pep),"\n",
        "keep.uniq: ",sum(keep.uniq),"\n")
  ## Applying filter
  if (!keepAll) {
    keep <- keep.na & keep.pep & keep.uniq & keep.int
  } else {
    keep <- TRUE
  }
  .featureData <- .featureData[keep,]
  .exprs <- .exprs[keep,]
  if (any(is.na(.featureData))) {
    whichCols <- apply(.featureData,2,function(x) any(is.na(x)))
    warning(paste("NA values in featureData column(s)",
                  names(.featureData)[whichCols]))
  }
  if (any(is.na(.exprs)))
    warning("NA values in .exprs.")
  ## Preparing return object slots and return new MSnSet
  .process <- new("MSnProcess",
                  processing=paste("Data loaded:",date(),
                    "from ispy result file."),
                  normalised=FALSE,
                  files=file)
  return(new("MSnSet",
             exprs=.exprs,
             featureData=new("AnnotatedDataFrame",data=.featureData),
             processingData=.process))
}
