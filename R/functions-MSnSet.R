
normalise.MSnSet <- function(object,method) {
  switch(method,
         max = div <- rowMax(exprs(object)),
         sum = div <- rowSums(exprs(object)))
  exprs(object) <- exprs(object)/div
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Normalised (",method,"): ",
                                              date(),
                                              sep=""))
  object@processingData@normalised <- TRUE
  return(object)
}

ratios.MSnSet <- function(object,log=FALSE) {
  r <- apply(exprs(object),1,getRatios,log)
  r <- t(r)
  rownames(r) <- featureNames(object)
  return(r)
}
