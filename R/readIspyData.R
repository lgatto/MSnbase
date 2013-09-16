readIspyData <- function(file = "ispy_results.tsv",
                         uniquePeps = TRUE,
                         pep = 0.05,
                         na.rm = TRUE,
                         min.int = 0,
                         reporters = 19:23,
                         ## skipFillUp = 24,
                         ## fillUp = FALSE,
                         keepAll = FALSE,
                         verbose = TRUE) {
  if (verbose)
    cat("Reading table\n")
  tab <- read.csv(file, header = TRUE, sep = "\t", fill = TRUE)
                  ## colClasses = c(
                  ##   "numeric",       # Index
                  ##   "factor",        # ProteinAccession 
                  ##   "character",     # ProteinDescription
                  ##   "numeric",       # ProteinQ.Value.uniquepeptidesonly.
                  ##   "numeric",       # ProteinType1Error.uniquepeptidesonly.
                  ##   "numeric",       # NumberOfUniquePeptides.FDR.1pc
                  ##   "character",     # PeptideSequence
                  ##   "numeric",       # PeptideParentProteins
                  ##   "character",     # MSFileID
                  ##   "numeric",       # PeptideType1ErrorForMSFile
                  ##   "factor",        # Fixed_Modifications
                  ##   "numeric",       # FixedModificationMassShift                    
                  ##   "factor",        # VariableModifications
                  ##   "numeric",       # VariableModificationMassShift
                  ##   "numeric",       # PrecursorMZ
                  ##   "character",     # MSScanID
                  ##   "character",     # DataBrowserLink ## was MascotQuery
                  ##   "numeric",       # PosteriorErrorProbability
                  ##   rep("numeric",length(reporters)),
                  ##   "numeric"))      # precursorrelativesignal
  names(tab) <- gsub("\\.\\.","pc", gsub("_","",names(tab)))
  names(tab)[ncol(tab)] <- "PrecursorRelativeSignal"
  ## This column was called MascotQuery, then changed to BrowserData
  ## - just ignore it
  ## tab$MascotQuery <- sub(",[0-0]+\\)&","",
  ##                        sub("=HYPERLINK\\(","",tab$MascotQuery))
  ## if (version==1)
  ##   tab$MascotQuery <- sub("^.+F0","F0",
  ##                          sub("\\.dat.+$",".dat",tab$MascotQuery))
  ## else
  ##   tab$DataBrowserLink <- sub("^.+F0","F0",
  ##                          sub("\\.dat.+$",".dat",tab$DataBrowserLink))
  .exprs <- as.matrix(tab[,reporters])
  .featureData <- tab[,-reporters] 
  ## if (fillUp) {
  ##   if (verbose)
  ##     cat("Filling up table\n")
  ##   if (verbose)
  ##     pb <- txtProgressBar(min = 0, max = ncol(.featureData), style = 3)
  ##   for (i in 1:ncol(.featureData)) {
  ##     .featureData[,i] <- fillUp(.featureData[,i])
  ##     if (verbose)
  ##       setTxtProgressBar(pb,i)
  ##   }
  ##   if (verbose)
  ##     close(pb)
  ## }
  if (verbose)
    cat("Filtering\n")
  ## Default filters keep all features
  keep.uniq <- keep.pep <- keep.na <- keep.int <- rep(TRUE,nrow(.exprs))
  ## Filtering on uniqueness of peptides
  if (uniquePeps) 
    keep.uniq <- as.numeric(.featureData$PeptideParentProteins) == 1
  ## Filtering on posterior error probability
  keep.pep <- as.numeric(.featureData$PosteriorErrorProbability) <= pep
  ## Remove features with NAs
  if (na.rm)
    keep.na <- apply(.exprs,1,function(x) !any(is.na(x)))
  rowsums <- rowSums(.exprs,na.rm=TRUE)
  keep.int <- rowsums >= min.int
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
  if (sum(keep) == 0)
    stop("No features left after filtering.")
  .exprs <- .exprs[keep,]
  .featureData <- .featureData[keep,]
  ## updating levels of factors
  .featureData$ProteinAccession <- .featureData$ProteinAccession[,drop=TRUE]
  .featureData$FixedModifications <- .featureData$FixedModifications[,drop=TRUE]
  .featureData$VariableModifications <- .featureData$VariableModifications[,drop=TRUE]
  if (any(is.na(.featureData))) {
    whichCols <- apply(.featureData,2,function(x) any(is.na(x)))
    nacolnames <- paste(names(.featureData)[whichCols], collapse = ", ")
    if (verbose)
      message(paste0("NA values in featureData column(s) ", nacolnames))
  }
  if (na.rm & any(is.na(.exprs)))
    warning("NA values in quantitation data.")
  ## Preparing return object slots and return new MSnSet
  .process <- new("MSnProcess",
                  processing = paste("Data loaded:",date(),
                    "from ispy result file."),
                  normalised = FALSE,
                  files = file)
  ans <- new("MSnSet",
             exprs=.exprs,
             featureData=new("AnnotatedDataFrame",data=.featureData))
  ans@processingData <- .process
  if (validObject(ans))
    return(ans)
}


readIspy15NData <- 
  readIspySilacData <-
  function(file = "ispy_results.tsv",
           uniquePeps = TRUE,
           pep = 0.05,
           na.rm = TRUE,
           min.int = 0,
           ...) {
    xx <- read.table(file, sep = "\t", header = TRUE, fill = TRUE, ...)
    ecols <- grep("precursor_intensity", names(xx))
    .exprs <- as.matrix(xx[, ecols])
    .fd <- new("AnnotatedDataFrame", data = xx[, -ecols])  
    ans <- new("MSnSet",
               exprs = .exprs,
               featureData = .fd)
    if (na.rm)
      ans <- filterNA(ans, pNA = 0)
    if (uniquePeps) 
      ans <- ans[featureData(ans)$Peptide_Parent_Proteins == 1,]  
    ans <- ans[featureData(ans)$Posterior_Error_Probability <= pep, ]
    ans <- ans[rowSums(exprs(ans), na.rm = TRUE) >= min.int, ]
    fData(ans) <- droplevels(fData(ans))
    if (validObject(ans))
      return(ans)
  }

