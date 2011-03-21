
normalise.MSnSet <- function(object,method) {
  switch(method,
         max = div <- rowMax(exprs(object)),
         sum = div <- rowSums(exprs(object)))
  exprs(object) <- exprs(object)/div
  object@process@processing <- c(object@process@processing,
                                 paste("Normalised (",method,"): ",
                                       date(),
                                       sep=""))
  object@process@normalised <- TRUE
  return(object)
}

ratios.MSnSet <- function(object,log=FALSE) {
  r <- apply(exprs(object),1,getRatios,log)
  r <- t(r)
  rownames(r) <- featureNames(object)
  return(r)
}
