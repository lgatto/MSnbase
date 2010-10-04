show.MSnQual <- function(object) {
  cat("Object of class \"",class(object),"\"\n",sep="")
  cat("Dimensions: ",dim(object@qc),"\n")
  show(object@metadata)
  show(object@process)
}
