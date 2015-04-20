msnsets <- function(object) object@x

MSnSetList <- function(x) .MSnSetList(x = x)

setMethod("show", "MSnSetList",
          function(object) {
              cat("Instance of class '", class(object), "' containig ",
                  length(object), " objects.\n", sep = "")
          })

setMethod("length", "MSnSetList", function(x) length(x@x))

setMethod("names", "MSnSetList", function(x) names(x@x))

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
          function(x, f, MARGIN = 1) {
              if (length(f) > 1)
                  stop("f must be of length 1.")
              if (MARGIN > 2 | MARGIN < 1)
                  stop("MARGIN must be 1 or 2.")
              if (MARGIN == 1) f <- factor(fData(x)[, f])
              else f <- factor(pData(x)[, f])
              split(x, f)
          })

setMethod("split", c("MSnSet", "factor"),
          function(x, f) {
              if (!length(f) %in% dim(x))
                  stop("length(f) not compatible with dim(x).")
              if (length(f) == nrow(x))
                  xl <- lapply(split(1:nrow(x), f), function(i) x[i, ])
              else ## length(f) == ncol(x)
                  xl <- lapply(split(1:ncol(x), f), function(i) x[, i])
              .MSnSetList(x = xl)
          })

## setMethod("unsplit", "MSnSetList",
##           function(value) {
##           })
