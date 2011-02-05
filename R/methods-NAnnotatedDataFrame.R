##################################################################
## Methods for NAnnotatedDataFrame class
setMethod("dim", "NAnnotatedDataFrame", function(x) {
  d <- c(dim(pData(x)),
         length(x@multiplex))
  names(d) <- c(dimLabels(x),"multipleNames")
  return(d)
})

setMethod("show",
          signature=signature(object="NAnnotatedDataFrame"),
          function(object) {
            callNextMethod(object)
            cat("  Multiplexing: ")
            if (length(object@multiplex)>0)
                cat(object@multiplex," - ",object@multiLabels,"\n",sep="")
            else
              cat("none\n")
          })


setMethod("multiplex", "NAnnotatedDataFrame",
          function(object) object@multiplex)

setMethod("multiLabels", "NAnnotatedDataFrame",
          function(object) object@multiLabels)
