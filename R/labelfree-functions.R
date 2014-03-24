# SIp Spectra Index per Peptide
# Griffin et al. 2010, PMID: 20010810
#
SI <- function(object, method=c("SIp", "SIgip", "SInp"), verbose=TRUE) {
  ## remove all missing identification data
  object <- removeNoId(object)

  ## remove all identification data with multiple assignments
  keep <- fData(object)$npsm == 1
  object <- object[keep, ]

  method <- match.arg(method)

  ## group by peptides seq
  groups <- as.factor(fData(object)$pepseq)
  ## SIp
  object <- combineFeatures(object, groupBy=groups, fun="sum", verbose=verbose)

  if (method %in% c("SIgip", "SInp")) {
    exprs(object) <- exprs(object)/colSums(exprs(object))
  }

  if (method == "SInp") {
    peplen <- fData(object)$len
    exprs(object) <- t(t(exprs(object))/peplen)
  }

  return(object)
}

