############################################################################
##                          DEPRECATED
## NAnnotatedDataFrame: As Biobase's AnnotatedDataFrame, it is composed of
## a data.frame, with annotations about columns named
## in the data slot contained in the metadata slot.
## In addition, it contains a multiplex slot to make explicite that
## the AnnotatedDataFrame is applied to a set of mulitplexed tags.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("NAnnotatedDataFrame",
         representation(multiplex = "numeric",
                        multiLabels = "character"),
         contains = c("AnnotatedDataFrame"),
         prototype = prototype(
              new("Versioned", versions = list(NAnnotatedDataFrame="0.0.3")),
             multiplex = 1,
             multiLabels = "Single run"),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@multiLabels) != object@multiplex)
                 msg <- validMsg(msg, "Number of multiplex does not match it's labels.")
             if (is.null(msg)) TRUE
             else msg
         })

setMethod("dim", "NAnnotatedDataFrame",
          function(x) {
              d <- c(dim(pData(x)),
                     x@multiplex)
              names(d) <- c(dimLabels(x), "multiNames")
              return(d)
          })

setMethod("show",
          signature = signature(object = "NAnnotatedDataFrame"),
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
