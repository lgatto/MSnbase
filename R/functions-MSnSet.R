#show.MSnSet <- function(object) {
#  cat("Object of class \"",class(object),"\"\n",sep="")
#  cat("Reporter ions:",names(object@quant),"\n")
#  cat(length(unique(object@precursors)),"precursors:",head(object@precursors,n=3),"...\n")
#  cat(length(object@precursors),"features:",head(object@features,n=3),"...\n")
#  show(object@metadata)
#  show(object@process)
#}

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

ratios.MSnSet <- function(object) {
  r <- apply(exprs(object),1,getRatios)
  r <- t(r)
  rownames(r) <- featureNames(object)
  return(r)
}
