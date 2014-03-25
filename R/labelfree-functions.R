# SIp Spectra Index per Peptide
# Griffin et al. 2010, PMID: 20010810
SI <- function(object, method=c("SIp", "SIgip", "SInp"), verbose=TRUE) {
  object <- .removeNoIdAndMultipleAssignments(object)

  method <- match.arg(method)

  ## group by protein
  groups <- as.factor(fData(object)$accession)
  ## SIp
  object <- combineFeatures(object, groupBy=groups, fun="sum", verbose=verbose)

  if (method %in% c("SIgip", "SInp")) {
    ## group by peptide 
    groups <- as.factor(fData(object)$pepseq)
    .exprs <- exprs(object)
    exprs(object) <- .exprs/ave(as.vector(.exprs), groups, FUN=sum)
  }

  if (method == "SInp") {
    ## divide by protein length
    protlen <- fData(object)$len
    exprs(object) <- t(t(exprs(object))/protlen)
  }

  return(object)
}

# (N)SAF (Normalised) Spectral Abundance Factor
# Paoletti et al. 2006, PMID: 17138671
SAF <- function(object, method=c("SAF", "NSAF"), verbose=TRUE) {
  object <- .removeNoIdAndMultipleAssignments(object)

  method <- match.arg(method)
  
  ## group by protein
  groups <- as.factor(fData(object)$accession)

  object <- combineFeatures(object, groupBy=groups, fun=length, verbose=verbose)

  ## divide by protein length
  protlen <- fData(object)$len
  exprs(object) <- t(t(exprs(object))/protlen)

  if (method == "NSAF") {
    ## normalize by all proteins of a complex (we use all values because we
    ## have no complex information)
    exprs(object) <- exprs(object)/colSums(exprs(object))
  }

  return(object)
}

.removeNoIdAndMultipleAssignments <- function(object) {
  ## remove all missing identification data
  object <- removeNoId(object)

  ## remove all identification data with multiple assignments
  keep <- fData(object)$npsm == 1
  object <- object[keep, ]
  
  return(object)
}
