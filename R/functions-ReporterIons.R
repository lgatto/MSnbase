show_ReporterIons <- function(object) {
  cat("Object of class \"", class(object), "\"\n",sep = "")
  cat(object@name, ": '", object@description, "' with ", length(object@mz),
      " reporter ions\n", sep = "")
  if (length(object@mz) > 0) {
      for (i in 1:length(object@mz))
          cat(" - [", object@reporterNames[i], "] ", object@mz[i], " +/- ",
            object@width, " (",object@col[i],")\n", sep = "")
  }
}


"[.ReporterIons" <- function(x,i) {
  msg <- validMsg(NULL, NULL)
  if (max(i)>length(x@mz) | min(i)<1)
    stop("subscript out of bonds")
  return(new("ReporterIons",
             name=paste(x@name,"[",paste(i,collapse=":"),"]",sep=""),
             description=paste("subset of",x@description,sep=" "),
             reporterNames=x@reporterNames[i],
             mz=x@mz[i],
             col=x@col[i],
             width=x@width))
}
