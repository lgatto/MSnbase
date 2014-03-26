##################################################################
## Methods for NAnnotatedDataFrame class
setMethod("dim", "NAnnotatedDataFrame", function(x) {
  d <- c(dim(pData(x)),
         x@multiplex)
  names(d) <- c(dimLabels(x), "multiNames")
  return(d)
})

setMethod("show",
          signature=signature(object="NAnnotatedDataFrame"),
          function(object) {
            callNextMethod(object)
            cat("  Multiplexing: ")
            if (length(object@multiplex) > 0)
                cat(object@multiplex, "-", object@multiLabels, "\n")
            else
              cat("none\n")
            invisible(NULL)
          })


setMethod("multiplex", "NAnnotatedDataFrame",
          function(object) object@multiplex)

setMethod("multiLabels", "NAnnotatedDataFrame",
          function(object) object@multiLabels)
