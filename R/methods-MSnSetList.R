msnsets <- function(object) object@x
objlog <- function(object) object@log

MSnSetList <-
    function(x, log = list(call = match.call()))
        .MSnSetList(x = x, log = log)

setMethod("show", "MSnSetList",
          function(object) {
              cat("Instance of class '", class(object), "' containig ",
                  length(object), " objects.\n", sep = "")
          })

setMethod("length", "MSnSetList", function(x) length(x@x))

setMethod("names", "MSnSetList", function(x) names(x@x))

setReplaceMethod("names", "MSnSetList",
          function(x, value) {
                     names(x@x) <- value
                     x
          })

setMethod("[", c("MSnSetList", "ANY", "missing", "missing"),
          function(x, i, j = "missing", drop = "missing")
              .MSnSetList(x = msnsets(x)[i]))

setMethod("[[", c("MSnSetList", "ANY", "missing"),
          function(x, i, j = "missing", drop = "missing") {
            if (length(i) != 1)
              stop("subscript out of bounds")
            msnsets(x)[[i]]
        })

setMethod("split", c("MSnSet", "character"),
          function(x, f) {
              if (length(f) != 1) stop("Character must be of lenght one.")
              if (f %in% varLabels(featureData(x))) {
                  f <- fData(x)[, f]
              } else if (f %in% varLabels(phenoData(x))) {
                  f <- pData(x)[, f]
              } else {
                  stop(f, " not found in any feature/phenodata variables.")
              }
              if (!is.factor(f)) f <- factor(f)
              split(x, f)
          })
          
setMethod("split", c("MSnSet", "factor"),
          function(x, f) {
              if (!length(f) %in% dim(x))
                  stop("length(f) not compatible with dim(x).")
              if (length(f) == nrow(x))
                  xl <- lapply(split(featureNames(x), f = f),
                               function(i) x[i, ])
              else ## length(f) == ncol(x)
                  xl <- lapply(split(sampleNames(x), f = f),
                               function(i) x[, i])
              MSnSetList(x = xl,
                         log = list(call = match.call(),
                             dims = dim(x),
                             f = f))
          })

setMethod("lapply", "MSnSetList",
          function(X, FUN, ...) {
              ans <- lapply(msnsets(X), FUN, ...)
              if (listOf(ans, "MSnSet"))
                  ans <- MSnSetList(ans)
              ans
          })

setMethod("unsplit", c("MSnSetList", "factor"),
          function(value, f) {
              len <- length(f)
              ## along what dimensions should we combine?
              ## (1) along rows
              dims1 <- c(ncol(value[[1L]]),
                         sum(unlist(lapply(value, nrow))))
              ## (2) along cols
              dims2 <- c(nrow(value[[1L]]),
                         sum(unlist(lapply(value, ncol))))
              if (!len %in% c(dims1, dims2))
                  stop(paste("length(f) is not compatible",
                             "with the object to be unsplit."))
              ans <- Reduce(combine, msnsets(value))
              if (len %in% dims1) {
                  xi <- lapply(value, featureNames)
                  xi <- unsplit(xi, f)
                  ans <- ans[xi, ]
              } else {
                  xi <- lapply(value, sampleNames)
                  xi <- unsplit(xi, f)
                  ans <- ans[, xi]
              }
              ans
          })

##
## un-exported utils
##

.sameNbCol <- function(x)
    length(unique(sapply(msnsets(x), ncol)) == 1)

.sameNbRow <- function(x)
    length(unique(sapply(msnsets(x), nrow)) == 1)
