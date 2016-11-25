setMethod("FeaturesOfInterest",
          c("character", "character", "missing"),
          function(fnames, description, ...) {
              .FeaturesOfInterest(fnames = fnames,
                                  description = description[1],
                                  date = date(),
                                  objpar = list())
          })

setMethod("FeaturesOfInterest",
          c("character", "character", "MSnSet"),
          function(fnames, description, object, ...) {
              fns <- featureNames(object)
              if (!all(fnames %in% fns)) {
                  fx <- !fnames %in% fns
                  stop(sum(fx), " feature(s) of interest absent from your object's feature names:\n   ", 
                       ifelse(sum(fx) > 5, 
                              paste(paste(head(fnames[fx], n = 5), collapse = ", "), ", ..."),
                              paste0(paste(fnames[fx], collapse = ", "), ".")))
              }
              .FeaturesOfInterest(fnames = fnames, 
                                  description = description[1], 
                                  date = date(),
                                  objpar = list(
                                      ncol = ncol(object),
                                      nrow = nrow(object),
                                      name = getVariableName(
                                          match.call(),
                                          "object"),
                                      digest = digest(object)))
              })

setMethod("foi", "FeaturesOfInterest",
          function(object, ...) object@fnames)

setMethod("description", "FeaturesOfInterest",
          function(object, ...) object@description)

setMethod("length", "FeaturesOfInterest",
         function(x) length(x@fnames))

setMethod("show", "FeaturesOfInterest",
          function(object) {
              if (length(object@objpar) > 0) cat("Traceable object")
              else cat("Object")
              cat(" of class \"", class(object), "\"\n", sep="")              
              cat(" Created on", object@date, "\n")
              cat(" Description:\n")
              cat(strwrap(object@description,
                          width = 0.65 * getOption("width"),
                          indent = 2, exdent = 2), sep = "\n")
              n <- length(object)
              cat(" ", n, " features of interest:\n", sep = "")
              if (n > 5) {
                  cat("  ", paste(object@fnames[1:2], collapse = ", "),
                      " ... ", paste(object@fnames[(n-1):n], collapse = ", "))
                  } else {
                      cat("  ", paste(object@fnames, collapse = ", "))
                  }
              cat("\n")
          })


setMethod("FoICollection", "missing",
          function(object, ...)
          .FoICollection(foic = list()))

setMethod("FoICollection", "list",
          function(object, ...) {
              val <- sapply(object, inherits, "FeaturesOfInterest")
              if (!all(val))
                  stop("The list must be composed of FeatureOfInterest instances")
              .FoICollection(foic = object)
          })

setMethod("length", "FoICollection",
         function(x) length(x@foic))

setMethod("lengths", "FoICollection",
          function(x, use.names = TRUE) {
              res <- sapply(x@foic, length)
              if (!use.names)
                  names(res) <- NULL
              res
          })

setMethod("show", "FoICollection",
          function(object)
          cat("A collection of ", length(object),
              " features of interest.\n", sep = ""))

setMethod("foi", "FoICollection",
          function(object, n) {
              if (missing(n)) {
                  return(object@foic)
              } else {
                  if (n > length(object))
                      stop("There are only ", length(object),
                           " available FeatureOfInterest instances.")
                  return(object@foic[n])
              }
          })


setMethod("addFeaturesOfInterest",
          c("FeaturesOfInterest", "FoICollection"),
          function(x, y) {
              .fois <- foi(y)
              check <- sapply(.fois, function(._foi) {
                  identical(sort(foi(._foi)), sort(foi(x)))
              })
              if (any(check)) {
                  message("The features of interest are already present.")
              } else {
                  y@foic <- c(y@foic, x)
              }
              return(y)
          })

setMethod("rmFeaturesOfInterest",
          c("FoICollection", "numeric"),
          function(object, i) {
              object@foic <- object@foic[-i]
              return(object)
          })

          
setMethod("description", "FoICollection",
          function(object, ...) sapply(foi(object), description))


## setMethod("fromIdentical",
##           c("FeaturesOfInterest", "FeaturesOfInterest"),
##           function(x, y) x@objpar$digest == y@objpar$digest)

## setMethod("fromIdentical",
##           c("FoICollection", "missing"),
##           function(x) {
##               dgsts <- sapply(foi(x), function(xx) xx@objpar$digest)
##               length(unique(dgsts)) == 1
##           })
              
## setMethod("fromEqual",
##           c("FeaturesOfInterest", "FeaturesOfInterest"),
##           function(x, y)
##           x@objpar$nrow == y@objpar$nrow & x@objpar$ncol == y@objpar$ncol)

## setMethod("fromEqual",
##           c("FoICollection", "missing"),
##           function(x) {
##               nc <- sapply(foi(x), function(xx) xx@objpar$ncol)
##               nr <- sapply(foi(x), function(xx) xx@objpar$nrow)
##               length(unique(nc)) == 1 & length(unique(nr)) == 1
##           })


## are any of the features in the FeaturesOfInterest instance present
## in the MSnSet (matrix) featureNames (rownames)?
setMethod("fnamesIn", c("FeaturesOfInterest", "MSnSet"),
          function(x, y, count = FALSE) fnamesIn(x, exprs(y), count))

setMethod("fnamesIn", c("FeaturesOfInterest", "data.frame"),
          function(x, y, count = FALSE) fnamesIn(x, as.matrix(y), count))
          
setMethod("fnamesIn", c("FeaturesOfInterest", "matrix"),
          function(x, y, count = FALSE) {             
              ans <- foi(x) %in% rownames(y)
              if (count) return(sum(ans))
              else return(any(ans))
          })


as.matrix.FoICollection <- function(x, ...) as(x, "matrix")

setAs("FoICollection", "matrix",
      function(from) {
          nms <- sapply(foi(from), description)
          names(nms) <- NULL          
          fns <- unique(unlist(lapply(foi(from), foi)))
          res <- matrix(0, ncol = length(nms), nrow = length(fns))
          rownames(res) <- fns
          colnames(res) <- nms
          for (k in foi(from))
              res[foi(k), description(k)] <- 1
          return(res)
      })

setMethod("names", "FoICollection", function(x) names(x@foic))

setReplaceMethod("names",
                 c("FoICollection", "character"),
                 function(x, value) {
                     names(x@foic) <- value
                     x
                 })
setMethod("[[", "FoICollection",
          function(x, i, j = "missing", drop = "missing")
              x@foic[[i]])

setMethod("[", "FoICollection",
          function(x, i, j = "missing", drop = "missing")
              FoICollection(x@foic[i]))
